# Assumptions and Limitations

This first version is a design-comparison model, not a clinical prediction
model. The executable baseline uses normalized carrier and cargo units because
the actual wet-lab dose unit, carrier size distribution, cargo loading, and
release yield have not been supplied.

## Assumptions used for execution

- IV bolus input is used: `A_blood(0) = D0`.
- The first reference system is a human exploratory context because the supplied
  project spreadsheet contains 70 kg human reference values.
- Missing kinetic constants are marked as exploratory assumptions in the
  parameter registry and are used only to exercise model structure.
- OCN does not alter tumor blood flow or EPR entry. It only changes local
  matrix binding, release, and OCN-bound terminal loss.
- YSA does not alter tumor blood flow or EPR entry. It only changes EphA2
  binding, unbinding, and receptor-mediated internalization.
- Non-intracellular terminal carrier loss removes unreleased payload into
  `A_cargo_lost`; only intracellular degradation produces active `A_cargo`.
- `alpha_ocn_cell_access` and the OCN-to-YSA transfer are phenomenological,
  unvalidated assumptions. Alpha is not a measured probability, does not create
  an `A_dual_bound` state, and does not alter tumor entry.
- Four designs share the same system parameters. No parameter is tuned to
  guarantee that OCN-YSA outperforms YSA-only.
- Forward binding and transfer totals are gross event counts; a carrier may be
  counted repeatedly after unbinding and rebinding.

## Limitations

- No supplied source gives measured CABOMB CaCO3 blood half-life.
- No supplied source gives CaCO3-specific liver or spleen uptake kinetics.
- No supplied source gives osteosarcoma-specific CABOMB tumor permeability or
  IFP.
- No supplied source gives OCN-hydroxyapatite association or dissociation
  kinetics after conjugation.
- No supplied source gives experimental YSA-EphA2 `Kd`, `kon`, `koff`, receptor
  density, or internalization rate.
- No supplied source gives pH-dependent CaCO3 degradation curves for the actual
  carrier.
- Docking scores are not converted into experimental binding constants.
- Paired Monte Carlo fractions depend on exploratory parameter distributions
  and are not probabilities of clinical success.

Chinese interpretation boundary: 这些比例是在探索性参数分布下得到的配对样本比例，
不是临床成功概率。

## V1.2 assumptions and limits

`Vmax_liver`, `Km_liver`, `k_liver_release`, `Vmax_spleen`, `Km_spleen`,
`k_spleen_release`, and `k_cross` are transparent, low-confidence exploratory
assumptions, not measured CABOMB values. Saturable liver/spleen uptake and
intact MPS release are unvalidated; MPS release is one candidate mechanism for
a slow terminal phase.

V1.2 is not a full PBPK model. It has no separate lung, kidney, normal-bone,
peripheral-tissue, arterial/venous, or transit compartments. Unmodeled organ
passage is not automatically terminal elimination, but is also not a
quantitative organ-distribution prediction. `C_blood` is normalized carrier per
apparent volume, not an absolute mass concentration. Sequential OCN-to-YSA
transfer is a model-structure assumption; fixed `R_EphA2_eff` is not a dynamic
receptor model. Model fractions are not probabilities of clinical success.

