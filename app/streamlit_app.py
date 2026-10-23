"""Streamlit interface for CABOMB V1.1, V1.2, and V1.3 exploratory modes."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from cabomb_model.designs import DESIGNS
from cabomb_model.fluxes import cumulative_flux_trajectories
from cabomb_model.metrics import summary_metrics, tumor_associated
from cabomb_model.model_modes import ModelModes, V1_1_MODES, V1_2_MODES
from cabomb_model.parameters import default_model_parameters, load_registry
from cabomb_model.scenarios import make_scenario
from cabomb_model.solver import simulate
from cabomb_model.v13_circulation import (
    SIZE_KEYS,
    v13_default_parameters,
    v13_detailed_metrics,
    v13_simulate,
    v13_summary_metrics,
    v13_tumor_associated,
)
from cabomb_model.v13_states import V13_STATES


st.set_page_config(page_title="CABOMB Model", layout="wide")

st.title("CABOMB PBPK-inspired ODE explorer")
st.warning(
    "Exploratory normalized model for design comparison only. "
    "Missing values are marked as assumptions and are not clinical predictions."
)

registry = pd.DataFrame(load_registry())

with st.sidebar:
    st.header("Scenario")
    model_version = st.selectbox(
        "Model mechanism version",
        ["V1.1 accessibility", "V1.2 sequential", "V1.3 explicit circulation"],
        help="V1.2 retains ten ODE states; V1.3 adds an explicit closed circulation (18 states).",
    )
    is_v13 = model_version == "V1.3 explicit circulation"
    if is_v13:
        params = v13_default_parameters()
        targeting_choice = st.selectbox(
            "Targeting mechanism",
            ["sequential_v1_2", "parallel_dual_access_v1_3"],
            index=0,
            help=(
                "parallel_dual_access_v1_3 adds a direct free-to-YSA route gated by "
                "alpha_ysa_direct_access, alongside the sequential free-to-OCN-to-YSA route."
            ),
        )
        modes = ModelModes("saturable_v1_2", targeting_choice, "explicit_v1_3")
    else:
        modes = V1_1_MODES if model_version == "V1.1 accessibility" else V1_2_MODES
        params = default_model_parameters()
        targeting_choice = modes.targeting_model

    st.caption(f"Active MPS model: {modes.mps_model}")
    st.caption(f"Active targeting model: {modes.targeting_model}")
    st.caption(f"Active circulation model: {modes.circulation_model}")
    design = st.selectbox("Carrier design", list(DESIGNS), index=list(DESIGNS).index("ocn-ysa"))
    epr = st.selectbox("EPR scenario", ["low", "medium", "high"], index=1)
    mps = st.selectbox("MPS clearance scenario", ["low", "medium", "high"], index=1)
    ifp = st.slider("Tumor IFP (mmHg)", 10.0, 45.0, float(params["IFP_tumor_project"]), 1.0)
    pcap_multiplier = st.slider("Permeability multiplier", 0.25, 3.0, 1.0, 0.05)
    if is_v13:
        size = st.selectbox("Size scenario", list(SIZE_KEYS), index=0)
    st.header("Local rates")
    params["k_ocn_on"] = st.slider("OCN association (h^-1)", 0.02, 0.8, float(params["k_ocn_on"]), 0.01)
    params["k_int"] = st.slider("YSA-bound internalization (h^-1)", 0.02, 1.0, float(params["k_int"]), 0.01)
    params["k_deg_int_acid"] = st.slider(
        "Acidic degradation (h^-1)",
        0.05,
        1.0,
        float(params["k_deg_int_acid"]),
        0.01,
    )
    if is_v13:
        st.header("V1.3 explicit-circulation rates")
        params["k_entry"] = st.slider(
            "k_entry (h^-1)", 0.005, 0.12, float(params["k_entry"]), 0.005
        )
        params["k_return"] = st.slider(
            "k_return (h^-1)", 0.02, 0.3, float(params["k_return"]), 0.01
        )
        params["Qtum"] = st.slider(
            "Tumor-feeding flow Qtum (L/min)", 0.01, 0.2, float(params["Qtum"]), 0.005
        )
        if targeting_choice == "parallel_dual_access_v1_3":
            params["alpha_ysa_direct_access"] = st.slider(
                "Direct free-to-YSA accessibility alpha_ysa_direct_access",
                0.0,
                1.0,
                float(params["alpha_ysa_direct_access"]),
                0.05,
                help=(
                    "Structural parallel dual-access assumption; 0 recovers the "
                    "sequential route, 1 gives full direct free-to-YSA access. Not a "
                    "measured binding constant."
                ),
            )
        st.info(
            "V1.3 circulation volumes and size factors are exploratory physiology "
            "references, not measured CABOMB values. OCN/YSA act only after "
            "extravasation and never change cardiac output or initial tumor supply."
        )
    elif modes == V1_1_MODES:
        params["k_ysa_on"] = st.slider(
            "YSA association (normalized receptor^-1 h^-1)",
            0.01,
            0.6,
            float(params["k_ysa_on"]),
            0.01,
        )
        st.header("Exploratory OCN accessibility")
        params["alpha_ocn_cell_access"] = st.slider(
            "OCN-bound cellular accessibility alpha",
            min_value=0.0,
            max_value=1.0,
            value=0.0,
            step=0.05,
            help=(
                "Phenomenological and unvalidated: effective YSA-EphA2 cellular access "
                "retained by OCN-bound carrier. This is not a measured probability, "
                "and it does not create a dual-bound state or change tumor entry."
            ),
        )
        if design != "ocn-ysa":
            st.info(
                "Alpha is inactive unless both OCN and YSA mechanism switches are enabled."
            )
    else:
        st.header("Exploratory V1.2 MPS and sequential targeting")
        params["Vmax_liver"] = st.slider(
            "Vmax_liver (normalized carrier / h)", 0.52, 6.24,
            float(params["Vmax_liver"]), 0.01,
        )
        params["Km_liver"] = st.slider(
            "Km_liver (normalized carrier / apparent L)", 0.2, 5.0,
            float(params["Km_liver"]), 0.01,
        )
        params["k_liver_release"] = st.slider(
            "k_liver_release (h^-1)", 0.0, 0.10,
            float(params["k_liver_release"]), 0.001,
        )
        params["Vmax_spleen"] = st.slider(
            "Vmax_spleen (normalized carrier / h)", 0.104, 2.60,
            float(params["Vmax_spleen"]), 0.01,
        )
        params["Km_spleen"] = st.slider(
            "Km_spleen (normalized carrier / apparent L)", 0.2, 5.0,
            float(params["Km_spleen"]), 0.01,
        )
        params["k_spleen_release"] = st.slider(
            "k_spleen_release (h^-1)", 0.0, 0.10,
            float(params["k_spleen_release"]), 0.001,
        )
        params["k_cross"] = st.slider(
            "k_cross (h^-1)", 0.0, 0.60,
            float(params["k_cross"]), 0.01,
        )
        st.info(
            "These saturable-MPS and sequential-targeting parameters are exploratory "
            "model assumptions. k_cross is not a measured probability or a validated "
            "receptor-binding constant."
        )

scenario = make_scenario(
    params,
    epr=epr,
    mps=mps,
    ifp_mmHg=ifp,
    pcap_tumor_cm_s=params["Pcap_tumor"] * pcap_multiplier,
)

if is_v13:
    result = v13_simulate(params, design=design, scenario=scenario, size=size, modes=modes)
    summary = v13_summary_metrics(result, params)

    top = st.columns(5)
    top[0].metric("Peak tumor-associated fraction", f"{summary['max_tumor_fraction']:.3f}")
    top[1].metric("Tumor AUC (normalized carrier x h)", f"{summary['tumor_auc']:.3f}")
    top[2].metric("Net extravasated fraction", f"{summary['net_extravasated_fraction']:.3f}")
    top[3].metric("Internalized (final fraction)", f"{summary['internalized_fraction_final']:.3e}")
    top[4].metric("Released cargo (final fraction)", f"{summary['released_cargo_fraction_final']:.3e}")
    st.metric("Maximum carrier mass-balance residual", f"{summary['max_abs_carrier_residual']:.1e}")
    st.caption(
        "V1.3 distinguishes venous, right-heart, pulmonary, left-heart, arterial, "
        "organ-branch, and tumor-vessel transit. Extravasation is a permeability-surface "
        "clearance on the tumor-vessel concentration."
    )

    circ = pd.DataFrame(
        [
            {"metric": "Fraction surviving first lung pass", "value": summary["fraction_surviving_first_lung_pass"]},
            {"metric": "Cumulative lung retention (fraction)", "value": summary["cumulative_lung_retention_fraction"]},
            {"metric": "Gross arterial delivery (fraction)", "value": summary["gross_arterial_delivery_fraction"]},
            {"metric": "Gross extravasated (fraction)", "value": summary["gross_extravasated_fraction"]},
            {"metric": "Net extravasated (fraction)", "value": summary["net_extravasated_fraction"]},
            {"metric": "Arterial AUC", "value": summary["arterial_auc"]},
            {"metric": "Tumor-vascular AUC", "value": summary["tumor_vascular_auc"]},
            {"metric": "Tumor-interstitial AUC", "value": summary["tumor_interstitial_auc"]},
        ]
    )
    st.subheader("V1.3 circulation and extravasation metrics")
    st.dataframe(circ, use_container_width=True)

    first_pass_result = v13_simulate(
        params, design=design, scenario=scenario, size=size, modes=modes, recirculate=False
    )
    detailed = v13_detailed_metrics(result, params, first_pass=first_pass_result)
    detail_rows = [
        {"metric": "First-pass lung arrival fraction", "value": detailed["first_pass_lung_arrival_fraction"]},
        {"metric": "First-pass lung retention fraction", "value": detailed["first_pass_lung_retention_fraction"]},
        {"metric": "First-pass systemic availability", "value": detailed["first_pass_systemic_availability"]},
        {"metric": "Tumor-vessel first-pass arrival fraction", "value": detailed["tumor_vascular_first_pass_arrival_fraction"]},
        {"metric": "Tumor first extravasation fraction", "value": detailed["tumor_first_extravasation_fraction"]},
        {"metric": "First-pass extravasation efficiency", "value": detailed["first_pass_extravasation_efficiency"]},
        {"metric": "Final tumor-associated stock / initial dose", "value": detailed["tumor_fraction_final"]},
        {"metric": "Apparent blood half-life (h)", "value": detailed["apparent_blood_half_life_h"]},
        {"metric": "OCN binding (gross)", "value": detailed["ocn_binding_fraction_gross"]},
        {"metric": "YSA binding (gross)", "value": detailed["ysa_binding_fraction_gross"]},
        {"metric": "OCN-to-YSA transfer (gross)", "value": detailed["ocn_to_ysa_transfer_fraction_gross"]},
        {"metric": "Total internalization / gross entry (specific + nonspecific)", "value": detailed["internalization_efficiency_gross"]},
        {"metric": "Specific internalization / gross entry", "value": detailed["ysa_specific_internalization_fraction_gross"]},
        {"metric": "Nonspecific internalization / gross entry", "value": detailed["nonspecific_internalization_fraction_gross"]},
        {"metric": "Cargo release efficiency (gross)", "value": detailed["cargo_release_efficiency_gross"]},
        {"metric": "Active cargo remaining (final)", "value": detailed["active_cargo_remaining_fraction_final"]},
    ]
    st.subheader("First-pass, delivery, and temporal metrics")
    st.dataframe(pd.DataFrame(detail_rows), use_container_width=True)

    trajectory = pd.DataFrame({"time_h": result.t})
    for state_name in V13_STATES:
        trajectory[state_name] = result.state(state_name)
    trajectory["tumor_associated_stock"] = v13_tumor_associated(result)

    st.subheader("Tumor-associated state stock")
    st.line_chart(trajectory.set_index("time_h")[["tumor_associated_stock"]], height=260)

    st.subheader("Explicit circulation trajectories")
    vascular_states = [
        "A_central_venous",
        "A_right_heart",
        "A_lung_vascular",
        "A_left_heart",
        "A_systemic_arterial",
        "A_liver_vascular",
        "A_spleen_vascular",
        "A_tumor_vascular",
        "A_rest_vascular",
    ]
    selected_states = st.multiselect(
        "States", V13_STATES, default=vascular_states + ["A_tumor_free", "A_internalized", "A_cargo"]
    )
    st.line_chart(trajectory.set_index("time_h")[selected_states], height=360)

    st.subheader("Four-design comparison (V1.3)")
    comparison_rows = []
    for design_name in DESIGNS:
        comparison_result = v13_simulate(params, design=design_name, scenario=scenario, size=size, modes=modes)
        comparison_rows.append(v13_summary_metrics(comparison_result, params))
    comparison = pd.DataFrame(comparison_rows)
    st.dataframe(
        comparison[
            [
                "design",
                "size",
                "max_tumor_fraction",
                "tumor_auc",
                "net_extravasated_fraction",
                "internalized_fraction_final",
                "released_cargo_fraction_final",
                "max_abs_carrier_residual",
            ]
        ],
        use_container_width=True,
    )
    st.warning(
        "V1.3 values are normalized exploratory amounts, not clinical predictions, "
        "and not a validated organ-distribution model."
    )
else:
    result = simulate(params, design=design, scenario=scenario, modes=modes)
    summary = summary_metrics(result, params)

    top = st.columns(5)
    top[0].metric(
        "Peak tumor-associated fraction",
        f"{summary['max_tumor_fraction']:.3f}",
    )
    top[1].metric(
        "Tumor AUC (normalized carrier x h)", f"{summary['tumor_auc']:.3f}"
    )
    top[2].metric(
        "Cumulative internalization (normalized carrier)",
        f"{summary['cumulative_internalization']:.4f}",
    )
    top[3].metric(
        "Cumulative cargo release (normalized cargo)",
        f"{summary['cumulative_cargo_release']:.4f}",
    )
    top[4].metric(
        "Cumulative OCN to YSA transfer (normalized carrier)",
        f"{summary['cumulative_ocn_to_ysa_transfer']:.4f}",
    )
    st.metric(
        "Maximum carrier mass-balance residual",
        f"{summary['max_abs_carrier_residual']:.1e}",
    )
    st.caption(
        "State stock is amount present at one time; exposure AUC integrates stock over "
        "time; cumulative flux integrates process throughput. Gross cumulative cargo "
        "release is not the same as active A_cargo remaining after loss."
    )
    if modes == V1_2_MODES:
        st.subheader("Cumulative V1.2 MPS processes")
        st.dataframe(
            pd.DataFrame(
                [
                    {
                        "organ": "liver",
                        "uptake": summary["cumulative_liver_mps_uptake"],
                        "intact_release_to_blood": summary["cumulative_liver_mps_release"],
                        "irreversible_degradation": summary["cumulative_liver_mps_degradation"],
                    },
                    {
                        "organ": "spleen",
                        "uptake": summary["cumulative_spleen_mps_uptake"],
                        "intact_release_to_blood": summary["cumulative_spleen_mps_release"],
                        "irreversible_degradation": summary["cumulative_spleen_mps_degradation"],
                    },
                ]
            ),
            use_container_width=True,
        )

    trajectory = pd.DataFrame({"time_h": result.t})
    for state_name in result.state_names:
        trajectory[state_name] = result.state(state_name)
    trajectory["tumor_associated_stock"] = tumor_associated(result)
    for name, values in cumulative_flux_trajectories(result, params).items():
        trajectory[name] = values

    st.subheader("Tumor-associated state stock")
    st.line_chart(trajectory.set_index("time_h")[["tumor_associated_stock"]], height=260)

    st.subheader("State trajectories")
    state_options = [
        "A_blood",
        "A_liver",
        "A_spleen",
        "A_tumor_free",
        "A_ocn_bound",
        "A_ysa_bound",
        "A_internalized",
        "A_elim_deg",
        "A_cargo",
    ]
    selected_states = st.multiselect("States", state_options, default=state_options[:7])
    st.line_chart(trajectory.set_index("time_h")[selected_states], height=360)

    st.subheader("Four-design comparison")
    rows = []
    for design_name in DESIGNS:
        comparison_result = simulate(
            params, design=design_name, scenario=scenario, modes=modes
        )
        rows.append(summary_metrics(comparison_result, params))
    comparison = pd.DataFrame(rows)
    st.dataframe(
        comparison[
            [
                "mps_model",
                "targeting_model",
                "design",
                "max_tumor_fraction",
                "tumor_auc",
                "internalized_fraction_final",
                "released_cargo_fraction_final",
                "cumulative_internalization",
                "cumulative_cargo_release",
                "cumulative_ocn_to_ysa_transfer",
                "cumulative_liver_mps_uptake",
                "cumulative_spleen_mps_uptake",
                "cumulative_liver_mps_release",
                "cumulative_spleen_mps_release",
                "cumulative_liver_mps_degradation",
                "cumulative_spleen_mps_degradation",
                "max_abs_carrier_residual",
            ]
        ],
        use_container_width=True,
    )
    st.warning(
        "Paired-sample fractions under exploratory parameter distributions are not probabilities of clinical success."
    )

st.subheader("Parameter registry")
st.dataframe(
    registry[
        [
            "parameter_symbol",
            "baseline_value",
            "lower_bound",
            "upper_bound",
            "unit",
            "value_type",
            "confidence_level",
            "notes",
        ]
    ],
    use_container_width=True,
    height=320,
)
