# CABOMB V1.3 Explicit-Circulation Model Specification

Status: specification for the optional V1.3 mode.

V1.3 is an **optional** mode. It does not modify, replace, or overwrite the
frozen V1.1/V1.2 ten-state model. The V1.1/V1.2 code path, its 128 tests, and
its frozen outputs remain byte-identical. V1.3 is selected with
`circulation_model = explicit_v1_3` (see below).

## 1. Purpose

V1.1/V1.2 use a single lumped `A_blood` compartment. V1.3 replaces it with an
explicit closed circulation so that the IV dose, first pulmonary pass, arterial
peak, organ-branch delivery, tumor-vascular arrival, transvascular
extravasation, interstitial retention, and venous return are distinguishable.
This addresses the "no explicit circulation" limitation recorded in
`docs/open_questions.md` without changing the tumor-targeting or MPS
mechanisms of V1.2.

## 2. Mode selection

The `ModelModes` dataclass gains a third field:

```python
ModelModes(mps_model, targeting_model, circulation_model="central_v1")
```

`CIRCULATION_MODELS = ("central_v1", "explicit_v1_3")`. The V1.3 baseline pair is

```python
V1_3_MODES = ModelModes("saturable_v1_2", "sequential_v1_2", "explicit_v1_3")
```

`central_v1` is the default and keeps V1.1/V1.2 behaviour unchanged. V1.3 keeps
the saturable MPS of V1.2; it only changes the circulation structure. V1.3 adds
a second, optional targeting mechanism, `parallel_dual_access_v1_3` (see §6),
which runs under the same explicit circulation.

`TARGETING_MODELS` therefore includes `"sequential_v1_2"` and
`"parallel_dual_access_v1_3"`.

## 3. State variables (18)

Circulation (vascular) states, unit = normalized carrier amount:

| Symbol | Definition |
| --- | --- |
| `A_central_venous` | Carrier in central/large systemic veins and vena cavae |
| `A_right_heart` | Carrier in right heart |
| `A_lung_vascular` | Carrier in pulmonary vascular bed |
| `A_left_heart` | Carrier in left heart |
| `A_systemic_arterial` | Carrier in systemic arteries (merged aortic/iliac/femoral/nutrient conduit, see §5) |
| `A_liver_vascular` | Carrier in liver sinusoidal blood |
| `A_spleen_vascular` | Carrier in splenic blood |
| `A_tumor_vascular` | Carrier in the osteosarcoma-feeding vessel (distal femoral) |
| `A_rest_vascular` | Carrier in all other tissue blood (muscle, brain, kidney, bone, GI, etc.) |

Tissue/organ states (kept from V1.2):

| Symbol | Definition |
| --- | --- |
| `A_liver` | Liver-MPS-associated carrier |
| `A_spleen` | Spleen-MPS-associated carrier |
| `A_tumor_free` | Free carrier in tumor interstitium |
| `A_ocn_bound` | OCN-retained carrier |
| `A_ysa_bound` | YSA-EphA2-bound carrier |
| `A_internalized` | Cellular-internalized carrier |
| `A_elim_deg` | Terminal eliminated/degraded carrier |
| `A_cargo` | Active released cargo |
| `A_cargo_lost` | Terminally removed cargo |

The carrier ledger consists of the 16 carrier-bearing states (all except
`A_cargo`, `A_cargo_lost`). The cargo ledger is `Y_cargo` times the 15 active
carrier states (all carrier states except `A_elim_deg`) plus `A_cargo` and
`A_cargo_lost`.

## 4. Assumptions (each is an explicit, low/medium-confidence assumption)

1. **Flow-limited transfer** `J = Q * C`, where `C_x = A_x / V_x`. No
   membrane-limited diffusion in the circulation; the only permeability-limited
   step is tumor-vessel-to-interstitium extravasation.
2. **Instantaneous IV bolus** enters `A_central_venous` at `t = 0` (initial
   condition `A_central_venous(0) = D0`), not a merged `A_blood`.
3. **No fake conduit states.** Aorta/iliac/femoral/nutrient arteries are not
   given independent ODE states with invented data; they are merged into
   `A_systemic_arterial`, whose outflow splits to the organ branches.
4. **Closed venous return.** Organ-vascular un-uptaken carrier returns to
   `A_central_venous`; nothing vanishes. Interstitial washout
   (`k_return * A_tumor_free`) returns to `A_central_venous` (lymphatic/venous
   drainage).
5. **Lung retention** is a terminal fraction `f_lung_ret` of pulmonary outflow
   (pulmonary-macrophage / transient sequestration), applied only at the lung
   compartment. It is size-dependent and below the >7 μm mechanical-filtration
   threshold for every CABOMB size. It never changes cardiac output.
6. **MPS uptake** is saturable and now sees the organ vascular concentration
   (`C_liver_vascular`, `C_spleen_vascular`) rather than a lumped blood
   concentration. Intact release returns to the organ vascular compartment.
7. **EPR/IFP/permeability** act only through `k_entry_eff` on the
   tumor-vascular-to-interstitium flux. They never alter cardiac output,
   pulmonary flow, or initial tumor-vascular delivery.
8. **OCN/YSA** never change blood flow, pulmonary transit, or initial tumor
   supply; OCN affects matrix retention and YSA affects post-arrival EphA2
   binding/internalization only (unchanged from V1.2 sequential structure).
9. **Lumped systemic clearance** `k_blood_clear` is applied to
   `A_rest_vascular` (the non-specialized tissue pool, including renal and
   other clearance), routing to `A_elim_deg`.

## 5. Flows and volumes

Branch flows (L/h) come from the project registry (R5): `Qco = 390`,
`Qliv = 87`, `Qspl = 11.4`, `Qtum = 2.7`, and
`Qrest = Qco - Qliv - Qspl - Qtum ≈ 288.9`. All circulations carry `Qco`
through the heart/lung; the systemic arterial outflow is the sum of the four
branch flows, which equals `Qco` by mass conservation.

Vascular volumes are listed in `docs/parameter_evidence_v1_3.md` §A.

## 6. Targeting (sequential V1.2 structure plus optional parallel dual access)

The baseline targeting structure is unchanged V1.2 sequential:
`A_tumor_free <-> A_ocn_bound -> A_ysa_bound -> A_internalized`. The only V1.3
circulation changes are (a) the entry term is a permeability-surface clearance
`PS_eff = k_entry_eff * size_extra * Vblood` acting on the tumor-vascular
concentration (`J_entry = PS_eff * C_tumor_vascular`), (b) `k_return` drains to
`A_central_venous`, and (c) internalization carries the `size_int` factor.
The four designs keep the same paths as V1.2 sequential mode.

The optional `parallel_dual_access_v1_3` mechanism adds a direct free-to-YSA
route alongside the sequential route, so that

```
A_tumor_free -> A_ocn_bound -> A_ysa_bound -> A_internalized   (sequential)
A_tumor_free --------------------> A_ysa_bound -> A_internalized   (direct)
```

The direct route is gated by `alpha_ysa_direct_access` (an assumed structural
accessibility parameter in `[0, 1]`; `0` recovers the sequential route exactly
and `1` gives full direct free-to-YSA access). Neither mechanism changes
cardiac output, pulmonary transit, EPR, permeability, or initial tumor-vessel
supply.

## 7. Size scenarios

A `size` selector (`100nm`, `200nm`, `500nm`, `800nm`, `2um`, `3um`) chooses
four size-dependent factors: lung retention `f_lung_ret`, MPS `size_mps`,
extravasation `size_extra`, and internalization `size_int` (see the evidence
table §B-E). These are the only size-dependent quantities.

## 8. Required V1.3 outputs

- Time curves for venous, RV, lung, LV, arterial, tumor-vascular, and
  interstitial (`A_tumor_free`) states, including the 0-5 min early pass and
  the 0-72 h horizon.
- Fraction surviving the first lung pass; fraction reaching tumor vascular;
  gross and net extravasated fraction.
- First-event (single-pass) delivery fractions, gross binding/internalization/
  release efficiencies, and temporal (peak, time-to-peak, 24/48/72 h, apparent
  blood half-life) metrics, all defined in `metric_dictionary.csv`.
- Final mass distribution across lung / liver / spleen / tumor / vascular /
  clearance.
- Arterial AUC, tumor-vascular AUC, tumor-interstitial AUC.
- Four-design ablation table and six-size comparison table.
- Sequential-vs-parallel mechanism comparison (both mechanisms, four designs,
  six sizes) and a feasibility-domain grid sweep of the targeting parameters.
- Local sensitivity, Monte Carlo (fixed seed), mass conservation,
  non-negativity, solver-stability, and fixed-seed reproducibility.
- Circulation + targeting mechanism figure (with the parallel direct route)
  and a mass-conservation figure.

## 9. Interpretation boundary

V1.3 remains a mechanistic exploration and design-comparison model. It is not
a clinical digital twin, efficacy predictor, human dose calculator, or
validated organ-distribution model. All V1.3 numbers are transparent
physiology references or exploratory assumptions (see the evidence table).
