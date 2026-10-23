# Equations

This file records the frozen first executable equation contract. The model uses
hours as the simulation time unit and carrier amount as the carrier state unit.
The first executable uses a normalized carrier dose for exploratory comparison
until a wet-lab dose unit is supplied.

## State variables

| Symbol | Definition | Unit |
| --- | --- | --- |
| `A_blood` | Circulating carrier in blood/central compartment | normalized carrier amount |
| `A_liver` | Liver-associated carrier | normalized carrier amount |
| `A_spleen` | Spleen-associated carrier | normalized carrier amount |
| `A_tumor_free` | Free carrier in tumor interstitium | normalized carrier amount |
| `A_ocn_bound` | OCN-retained carrier bound to hydroxyapatite/bone matrix | normalized carrier amount |
| `A_ysa_bound` | YSA-EphA2-bound carrier | normalized carrier amount |
| `A_internalized` | Cellular internalized carrier | normalized carrier amount |
| `A_elim_deg` | Eliminated or terminally degraded carrier | normalized carrier amount |
| `A_cargo` | Active released cargo still represented as available | normalized cargo amount |
| `A_cargo_lost` | Released, exported, degraded, or terminally removed cargo | normalized cargo amount |

## Tumor entry function

The ODE uses `k_entry_eff` in `h^-1`. Source values such as permeability, IFP,
EPR pore size, and blood flow are not used directly as first-order rates.

```text
permeability_factor = Pcap_tumor / Pcap_tumor_ref
ifp_factor = (1 + ifp_entry_sensitivity * IFP_ref)
             / (1 + ifp_entry_sensitivity * IFP_tumor)
k_entry_eff = k_entry * EPR_factor * permeability_factor * ifp_factor
```

Higher IFP therefore decreases entry all else equal, while higher effective
permeability or EPR scenario factor increases entry.

## ODE system

```text
dA_blood/dt =
  - (k_liver_uptake_eff + k_spleen_uptake_eff + k_entry_eff + k_blood_clear) A_blood
  + k_return A_tumor_free

dA_liver/dt =
  k_liver_uptake_eff A_blood - k_liver_elim A_liver

dA_spleen/dt =
  k_spleen_uptake_eff A_blood - k_spleen_elim A_spleen

dA_tumor_free/dt =
  k_entry_eff A_blood
  - k_return A_tumor_free
  - I_OCN k_ocn_on A_tumor_free + I_OCN k_ocn_off A_ocn_bound
  - I_YSA k_ysa_on R_EphA2_eff A_tumor_free + I_YSA k_ysa_off A_ysa_bound
  - k_int_ns A_tumor_free
  - k_tumor_deg A_tumor_free

dA_ocn_bound/dt =
  I_OCN k_ocn_on A_tumor_free
  - I_OCN k_ocn_off A_ocn_bound
  - I_OCN k_ocn_loss A_ocn_bound

dA_ysa_bound/dt =
  I_YSA k_ysa_on R_EphA2_eff A_tumor_free
  - I_YSA k_ysa_off A_ysa_bound
  - I_YSA k_int A_ysa_bound

dA_internalized/dt =
  I_YSA k_int A_ysa_bound
  + k_int_ns A_tumor_free
  - k_deg_int(pH) A_internalized

dA_elim_deg/dt =
  k_blood_clear A_blood
  + k_liver_elim A_liver
  + k_spleen_elim A_spleen
  + k_tumor_deg A_tumor_free
  + I_OCN k_ocn_loss A_ocn_bound
  + k_deg_int(pH) A_internalized

dA_cargo/dt =
  Y_cargo k_deg_int(pH) A_internalized
  - k_cargo_loss A_cargo

dA_cargo_lost/dt =
  Y_cargo (
    k_blood_clear A_blood
    + k_liver_elim A_liver
    + k_spleen_elim A_spleen
    + k_tumor_deg A_tumor_free
    + I_OCN k_ocn_loss A_ocn_bound
  )
  + k_cargo_loss A_cargo
```

The carrier ledger is separate from the cargo ledger:

```text
carrier_total =
  A_blood + A_liver + A_spleen + A_tumor_free + A_ocn_bound
  + A_ysa_bound + A_internalized + A_elim_deg

cargo_total =
  Y_cargo * (A_blood + A_liver + A_spleen + A_tumor_free + A_ocn_bound
             + A_ysa_bound + A_internalized)
  + A_cargo + A_cargo_lost
```

The first ledger must equal `D0` within numerical tolerance. The second must
equal `Y_cargo * D0` when cargo accounting is enabled.

## V1.1 OCN-to-YSA transfer

```text
J_ocn_to_ysa =
  I_OCN * I_YSA * alpha_ocn_cell_access
  * k_ysa_on * R_EphA2_eff * A_ocn_bound

dA_ocn_bound/dt = V1 OCN-bound terms - J_ocn_to_ysa
dA_ysa_bound/dt = V1 YSA-bound terms + J_ocn_to_ysa
```

The equal and opposite terms conserve carrier amount and do not represent an
explicit or hidden dual-bound state. `J_ocn_to_ysa` does not occur in the tumor
entry expression.

For each non-negative instantaneous process flux `J_q(t)`, cumulative throughput
is postprocessed on the saved time grid using composite trapezoidal integration:

```text
C_q(t_0) = 0
C_q(t_n) = integral from t_0 to t_n of J_q(t) dt
```

Required totals include tumor entry and return, OCN binding, YSA binding from
free carrier, OCN-to-YSA transfer, total and component internalization,
intracellular degradation, and cargo release. These cumulative process totals
are distinct from state stock and exposure AUC. In particular,
`cumulative_cargo_release` is gross historical release before cargo loss,
whereas `A_cargo(t)` is active released cargo remaining at time `t`.

## V1.2 optional mode equations

`C_blood = A_blood / Vblood` is normalized carrier per apparent volume, not an
absolute mass concentration. In the saturable mode:

```text
J_liver_MPS = mps_factor * Vmax_liver * C_blood / (Km_liver + C_blood)
J_spleen_MPS = mps_factor * Vmax_spleen * C_blood / (Km_spleen + C_blood)
J_liver_release = k_liver_release * A_liver
J_spleen_release = k_spleen_release * A_spleen
```

Release is an equal return of intact carrier to blood; organ degradation remains
an irreversible contribution to `A_elim_deg`. Saturable uptake and intact MPS
release are unvalidated CABOMB assumptions, and MPS release is one candidate
mechanism for a slow terminal phase.

### Saturable MPS right-hand-side replacements

When `mps_model = saturable_v1_2`, the MPS-dependent state equations are:

```text
dA_blood/dt =
  - J_liver_MPS - J_spleen_MPS
  + J_liver_release + J_spleen_release
  - k_entry_eff A_blood
  - k_blood_clear A_blood
  + k_return A_tumor_free

dA_liver/dt =
  J_liver_MPS - J_liver_release - k_liver_elim A_liver

dA_spleen/dt =
  J_spleen_MPS - J_spleen_release - k_spleen_elim A_spleen

dA_elim_deg/dt =
  k_blood_clear A_blood
  + k_liver_elim A_liver + k_spleen_elim A_spleen
  + k_tumor_deg A_tumor_free
  + I_OCN k_ocn_loss A_ocn_bound
  + k_deg_int(pH) A_internalized

dA_cargo_lost/dt =
  Y_cargo (
    k_blood_clear A_blood
    + k_liver_elim A_liver + k_spleen_elim A_spleen
    + k_tumor_deg A_tumor_free
    + I_OCN k_ocn_loss A_ocn_bound
  )
  + k_cargo_loss A_cargo
```

Thus each intact release flux is added to blood and subtracted from its organ
with equal magnitude. `k_liver_elim A_liver` and
`k_spleen_elim A_spleen` remain irreversible degradation fluxes into the
terminal carrier state, with the corresponding unreleased payload entered in
the lost-cargo ledger. J_liver_release and J_spleen_release do not enter
`A_elim_deg` or `A_cargo_lost`.

When `mps_model = linear_v1`, the original first-order uptake equations above
remain active and both intact-release fluxes are zero.

For strict sequential OCN-YSA targeting:

```text
J_ysa_from_free = I_YSA * (1 - I_OCN) * k_ysa_on * R_EphA2_eff * A_tumor_free
J_ocn_to_ysa = I_OCN * I_YSA * k_cross * R_EphA2_eff * A_ocn_bound

dA_tumor_free/dt =
  k_entry_eff A_blood
  - k_return A_tumor_free
  - I_OCN k_ocn_on A_tumor_free
  + I_OCN k_ocn_off A_ocn_bound
  - J_ysa_from_free
  + I_YSA k_ysa_off A_ysa_bound
  - k_int_ns A_tumor_free
  - k_tumor_deg A_tumor_free

dA_ocn_bound/dt =
  I_OCN k_ocn_on A_tumor_free
  - I_OCN k_ocn_off A_ocn_bound
  - I_OCN k_ocn_loss A_ocn_bound
  - J_ocn_to_ysa

dA_ysa_bound/dt =
  J_ysa_from_free
  + J_ocn_to_ysa
  - I_YSA k_ysa_off A_ysa_bound
  - I_YSA k_int A_ysa_bound
```

The transfer conserves carrier between `A_ocn_bound` and `A_ysa_bound`, does not
alter tumor entry, and is a sequential structure assumption. Fixed
`R_EphA2_eff` means this is not a dynamic receptor or receptor-saturation model.
For the OCN-YSA design, `J_ysa_from_free = 0`; its targeted forward path must
therefore pass through OCN retention before YSA-EphA2 binding. YSA dissociation
still returns carrier from `A_ysa_bound` to `A_tumor_free`.

The four designs follow these targeting paths in `sequential_v1_2`:

| Design | Active targeting path |
| --- | --- |
| unmodified | No OCN or YSA targeting path; tumor-free carrier can still return, degrade, or internalize nonspecifically. |
| OCN only | `A_tumor_free <-> A_ocn_bound`, with OCN loss; there is no YSA-bound path. |
| YSA only | `A_tumor_free <-> A_ysa_bound -> A_internalized` through direct free-to-YSA binding. |
| OCN-YSA | `A_tumor_free <-> A_ocn_bound -> A_ysa_bound -> A_internalized`; direct free-to-YSA binding is zero. |

## 2026-09-17 V1.3 reporting clarification

The V1.3 ODE is unchanged. Total cumulative internalization is the time integral
of `internalization + nonspecific_internalization`; its gross efficiency divides
this sum by cumulative `entry`, including repeat entry. Dose-normalized cumulative
specific, nonspecific and total internalization are reported separately from
`A_internalized(final) / D0`, which is the surviving internalized carrier stock.
Intracellular cargo release is `Y_cargo * integral(deg_int dt)`; final active cargo
also reflects `k_cargo_loss`. Neither quantity is calcium mass or free Ca2+.

The legacy `tumor_retention_efficiency` actually divides final tumor-associated
stock by gross cumulative entry in the same recirculating run, not first entry.
It is excluded from core reporting; the implementation remains for compatibility.
Report final tumor stock and cumulative entry separately. Tumor AUC integrates
carrier amount (unit: normalized carrier amount * h); it is not divided by D0.
