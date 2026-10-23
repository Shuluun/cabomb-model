# CABOMB V1.3 Explicit-Circulation Equations

This is the equation contract for the optional `explicit_v1_3` circulation
mode. Time unit is hours; carrier amount is the state unit; the carrier dose is
normalized. V1.3 uses `mps_model = saturable_v1_2` and
`targeting_model = sequential_v1_2`, so its tumor-targeting and saturable-MPS
terms are the V1.2 forms re-expressed over the circulation states.

## Vascular concentrations

```text
C_x = A_x / V_x      for each vascular compartment x
```

## Circulation fluxes (flow-limited, J = Q * C)

```text
Qrest = Qco - Qliv - Qspl - Qtum          (mass-conservation remainder)

J_cv_to_rh   = Qco * C_cv
J_rh_to_lung = Qco * C_rh
J_lung_pass  = Qco * C_lung * (1 - f_lung_ret)
J_lung_ret   = Qco * C_lung * f_lung_ret        (terminal)
J_lh_to_sa   = Qco * C_lh

J_sa_to_liver  = Qliv * C_sa
J_sa_to_spleen = Qspl * C_sa
J_sa_to_tumor  = Qtum * C_sa
J_sa_to_rest   = Qrest * C_sa

J_liver_venous = Qliv * C_liver_vascular
J_spleen_venous = Qspl * C_spleen_vascular
J_tumor_venous = Qtum * C_tumor_vascular
J_rest_venous  = Qrest * C_rest_vascular
```

## MPS and tumor-entry fluxes

```text
J_liver_MPS = mps_factor * size_mps * Vmax_liver * C_liver_vascular
              / (Km_liver + C_liver_vascular)
J_spleen_MPS = mps_factor * size_mps * Vmax_spleen * C_spleen_vascular
               / (Km_spleen + C_spleen_vascular)
J_liver_release = k_liver_release * A_liver
J_spleen_release = k_spleen_release * A_spleen

PS_eff = k_entry_eff * size_extra * Vblood       (permeability-surface clearance)
J_entry   = PS_eff * C_tumor_vascular = PS_eff * A_tumor_vascular / V_tumor_vascular
J_return  = k_return * A_tumor_free
J_blood_clear = k_blood_clear * A_rest_vascular        (terminal)
```

`k_entry_eff` is computed by the existing `compute_k_entry` function
(EPR factor * permeability factor * IFP factor). EPR/IFP/permeability appear
nowhere else. `J_entry` is a permeability-surface clearance `PS_eff` acting on
the tumor-vascular concentration, so it is flow- and permeability-limited with
extraction fraction `E = PS_eff / (Qtum + PS_eff) ∈ [0, 1)`, and it reproduces
the V1.2 lumped tumor-entry scale in the well-mixed limit.

## Targeting fluxes (sequential_v1_2 and parallel_dual_access_v1_3)

```text
J_ocn_on         = I_OCN * k_ocn_on * A_tumor_free
J_ocn_off        = I_OCN * k_ocn_off * A_ocn_bound
J_ocn_loss       = I_OCN * k_ocn_loss * A_ocn_bound
J_ocn_to_ysa     = I_OCN * I_YSA * k_cross * R_EphA2_eff * A_ocn_bound
J_ysa_off        = I_YSA * k_ysa_off * A_ysa_bound
J_int            = I_YSA * k_int * size_int * A_ysa_bound
J_int_ns         = k_int_ns * size_int * A_tumor_free
J_deg_int        = k_deg_int(pH) * A_internalized
J_tumor_deg      = k_tumor_deg * A_tumor_free
```

The free-to-YSA flux selects between the two mechanisms:

```text
sequential_v1_2:
    J_ysa_from_free = I_YSA * (1 - I_OCN) * k_ysa_on * R_EphA2_eff * A_tumor_free

parallel_dual_access_v1_3:
    J_ysa_from_free = I_YSA * [(1 - I_OCN) + I_OCN * alpha_ysa_direct_access]
                      * k_ysa_on * R_EphA2_eff * A_tumor_free
```

`alpha_ysa_direct_access ∈ [0, 1]` is the direct free-to-YSA accessibility of
OCN-bearing carrier. `alpha_ysa_direct_access = 0` recovers the sequential route
exactly; `1` gives full direct free-to-YSA access in parallel with the
sequential free-to-OCN-to-YSA route. Neither mechanism changes any circulation,
entry, or extravasation flux.

## ODE system

```text
dA_cv/dt  = -J_cv_to_rh + J_liver_venous + J_spleen_venous
            + J_tumor_venous + J_rest_venous + J_return

dA_rh/dt  = J_cv_to_rh - J_rh_to_lung

dA_lung/dt = J_rh_to_lung - J_lung_pass - J_lung_ret

dA_lh/dt  = J_lung_pass - J_lh_to_sa

dA_sa/dt  = J_lh_to_sa - J_sa_to_liver - J_sa_to_spleen
            - J_sa_to_tumor - J_sa_to_rest

dA_liv_v/dt = J_sa_to_liver + J_liver_release - J_liver_venous - J_liver_MPS

dA_spl_v/dt = J_sa_to_spleen + J_spleen_release - J_spleen_venous - J_spleen_MPS

dA_tum_v/dt = J_sa_to_tumor - J_tumor_venous - J_entry

dA_rest_v/dt = J_sa_to_rest - J_rest_venous - J_blood_clear

dA_liver/dt = J_liver_MPS - J_liver_release - k_liver_elim * A_liver

dA_spleen/dt = J_spleen_MPS - J_spleen_release - k_spleen_elim * A_spleen

dA_tumor_free/dt = J_entry - J_return - J_ocn_on + J_ocn_off
                   - J_ysa_from_free + J_ysa_off - J_int_ns - J_tumor_deg

dA_ocn_bound/dt = J_ocn_on - J_ocn_off - J_ocn_loss - J_ocn_to_ysa

dA_ysa_bound/dt = J_ysa_from_free + J_ocn_to_ysa - J_ysa_off - J_int

dA_internalized/dt = J_int + J_int_ns - J_deg_int

dA_elim_deg/dt = J_lung_ret + J_blood_clear
                 + k_liver_elim * A_liver + k_spleen_elim * A_spleen
                 + J_tumor_deg + J_ocn_loss + J_deg_int

dA_cargo/dt = Y_cargo * J_deg_int - k_cargo_loss * A_cargo

dA_cargo_lost/dt = Y_cargo * (J_lung_ret + J_blood_clear
                              + k_liver_elim * A_liver + k_spleen_elim * A_spleen
                              + J_tumor_deg + J_ocn_loss)
                   + k_cargo_loss * A_cargo
```

## Conservation ledgers

```text
carrier_total = A_cv + A_rh + A_lung + A_lh + A_sa + A_liv_v + A_spl_v
                + A_tum_v + A_rest_v + A_liver + A_spleen + A_tumor_free
                + A_ocn_bound + A_ysa_bound + A_internalized + A_elim_deg
carrier_total = D0          (within numerical tolerance)

active_carrier = carrier_total - A_elim_deg
cargo_total = Y_cargo * active_carrier + A_cargo + A_cargo_lost
cargo_total = Y_cargo * D0  (within numerical tolerance)
```

Every flow-limited circulation flux appears as an equal and opposite pair, so
the carrier ledger closes. The cargo ledger closes because every terminal
carrier flux carries `Y_cargo` payload into `A_cargo_lost`, and only
intracellular degradation releases payload into `A_cargo`.

## Initial condition and time grid

```text
A_cv(0) = D0; all other states = 0
```

The solver uses LSODA with `rtol = 1e-8`, `atol = 1e-10` and a dense early
time grid: `linspace(0, 0.05 h, 400)` (≈ 0.45 s spacing to resolve the first
pulmonary pass and arterial peak) followed by `linspace(0.05, 72 h, 240)[1:]`.

## First-pass and extravasation metrics

```text
f_surviving_first_lung_pass = 1 - f_lung_ret
cumulative_lung_ret = ∫ J_lung_ret dt
cumulative_arterial_delivery = ∫ J_sa_to_tumor dt        (gross, can exceed D0 via recirculation)
cumulative_entry = ∫ J_entry dt                          (gross extravasation)
cumulative_return = ∫ J_return dt
net_extravasated = cumulative_entry - cumulative_return
f_truly_extravasated = net_extravasated / D0
```

All cumulative integrals are composite-trapezoid postprocessing on the saved
time grid, matching the V1.2 convention.

### First-event (single-pass) fractions vs gross throughput

A single-pass run (`recirculate = False`) diverts venous return to a sink so no
particle recirculates. Its cumulative fluxes are therefore unique-particle
**first-event** fractions, never exceeding `D0`:

```text
first_pass_lung_arrival_fraction   = cumulative_rh_to_lung / D0
first_pass_lung_retention_fraction = cumulative_lung_ret / D0
first_pass_systemic_availability   = cumulative_lh_to_sa / D0
tumor_vascular_first_pass_arrival  = cumulative_sa_to_tumor / D0
tumor_first_extravasation_fraction = cumulative_entry / D0
first_pass_extravasation_efficiency = cumulative_entry / cumulative_sa_to_tumor
```

By contrast, the recirculating-run cumulative fluxes are **gross** throughput
and may exceed `D0` or 100 %; binding/internalization/release efficiencies
normalized to gross extravasation are process efficiencies, not binding
probabilities. The full metric dictionary (numerator, denominator, unit, and
stock/gross/first-event classification) is in
`outputs/v1_3/tables/metric_dictionary.csv`.
