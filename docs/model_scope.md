# Model Scope for CABOMB V1.1

This first version is a PBPK-inspired compartmental ODE specification for design
comparison and mechanistic exploration. It is not a clinical dose-prediction
model, not a patient digital twin, and not a 3D vascular simulation.

V1.1 retains exactly ten ODE states. It does not add `A_dual_bound`, cumulative
counter states, a complete PBPK network, lung or normal-bone compartments,
spatial diffusion, or a killing model. Four carrier designs share the same
system parameters. OCN and YSA remain local post-entry mechanisms and do not
alter blood flow or tumor entry.

V1.1 separately reports state stock (amount present at a time), exposure AUC
(stock integrated over time), and cumulative flux (gross process throughput
integrated over time). The exploratory OCN-to-YSA transfer is phenomenological
and unvalidated; alpha zero preserves accepted V1 behavior.

## Included in Version 1

1. Intravenous systemic circulation of CaCO3 carriers.
   - Basis: project PPT and `igem_model_cases.pdf` frame the primary question as
     how many microspheres reach the osteosarcoma region after IV injection.
   - Practical reason: a central circulation state is the smallest useful entry
     point for mass balance and dose tracking.

2. Liver and spleen MPS/RES uptake.
   - Basis: `nanoparticle_PBPK.pdf`, `delivering_nanomedicine.pdf`,
     `modeling_guidance.docx`, and the project spreadsheet all identify
     macrophage/RES/MPS uptake, especially liver and spleen, as a major
     nanoparticle clearance route.
   - Practical reason: ignoring this would overstate circulating carrier and
     tumor availability.

3. Tumor extravasation as an uncertain EPR/IFP/permeability process.
   - Basis: `EPR_review_2022.pdf` emphasizes EPR heterogeneity and limited
     clinical predictability; `tumor_delivery_2011.pdf` and
     `delivering_nanomedicine.pdf` identify permeability, pore size, IFP, ECM,
     and lymphatic drainage as barriers.
   - Practical reason: the first model should compare low, medium, and high EPR
     scenarios rather than assume guaranteed tumor entry.

4. Tumor interstitial free carrier and local retention.
   - Basis: tumor delivery reviews distinguish extravasation, interstitial
     transport, binding barriers, and retention.
   - Practical reason: this separates "arrived near tumor" from "retained where
     it can interact locally".

5. OCN-mediated matrix binding/retention.
   - Basis: project PPT describes OCN affinity for hydroxyapatite/bone matrix.
   - Practical reason: OCN should change local retention after tumor-region
     arrival, not systemic blood flow.

6. YSA-EphA2 binding and receptor-mediated internalization.
   - Basis: project PPT describes YSA as recognizing EphA2 on osteosarcoma cells;
     tumor delivery reviews describe ligand-mediated binding and endocytosis as
     post-arrival processes.
   - Practical reason: this lets OCN-only, YSA-only, and OCN-YSA designs be
     compared mechanistically.

7. Acid-dependent CaCO3 degradation and cargo release.
   - Basis: project PPT and `igem_model_cases.pdf` describe acidic
     endosome/lysosome degradation and release of payload; reviews discuss acidic
     tumor/endosomal environments and intracellular trafficking.
   - Practical reason: release must be downstream of internalization/degradation,
     not identical to tumor arrival.

8. Explicit degraded/eliminated carrier state.
   - Basis: AGENTS.md requires mass-balance checks.
   - Practical reason: eliminated/degraded carrier must be tracked separately so
     total carrier dose can be audited.

## Excluded from Version 1

1. Full 3D vascular geometry, CFD, Navier-Stokes, Langevin particle tracking, and
   margination modeling.
   - Basis: PPT slide 5 image identifies these as micro-scale methods, but the
     mission specifies a PBPK-inspired ODE first version.
   - Practical reason: these require vessel geometry, local flow fields, and
     particle-wall interaction data not supplied here.

2. Active targeting as a force, navigation mechanism, or direct increase in
   tumor blood flow.
   - Basis: AGENTS.md and project PDF/PPT text explicitly reject this. Active
     targeting enters local binding, retention, and internalization only.
   - Practical reason: otherwise the model would confound tumor entry with
     post-arrival cell/matrix interactions.

3. Full immune-cell mechanistic model.
   - Basis: project PDF says to omit detailed immune-cell phagocytosis initially;
     literature supports MPS importance but does not provide CaCO3-specific rates.
   - Practical reason: liver/spleen uptake can be represented by effective rates
     until experiments or literature provide detailed kinetics.

4. Full MNAzyme sequence-design platform and intracellular killing model.
   - Basis: project PDF lists MNAzyme design as a separate later workstream.
   - Practical reason: the first delivery model can track released cargo, but
     cannot predict sequence-specific MNAzyme efficacy without target/trigger
     sequences and validation data.

5. Docking-score-based experimental `Kd` prediction.
   - Basis: AGENTS.md forbids silently converting docking scores into
     experimental `Kd`; project PDF mentions this idea but it remains unvalidated.
   - Practical reason: docking can guide qualitative ranking, but cannot serve as
     a quantitative binding parameter without calibration.

6. Human clinical prediction or patient-specific dosing.
   - Basis: EPR and PBPK literature emphasize species, tumor, patient, and
     material heterogeneity.
   - Practical reason: supplied data are insufficient for validated human
     prediction.

## Questions Version 1 Can Answer

- Under the same dose and scenario assumptions, which design retains more carrier
  in tumor-associated states: unmodified, OCN-only, YSA-only, or OCN-YSA?
- How sensitive is qualitative design ranking to MPS clearance and EPR strength?
- Which missing parameters most affect tumor free, matrix-bound, receptor-bound,
  internalized, degraded, and released-cargo outputs?
- What experimental measurements would most improve the next model iteration?

## Questions Version 1 Cannot Answer

- A clinical human dose, efficacy, safety margin, or tumor cure probability.
- Patient-specific or tumor-location-specific delivery without imaging or
  patient data.
- True 3D carrier trajectories, margination, collision with vessel walls, or
  local extravasation hotspots.
- Experimental `Kd`, `kon`, `koff`, internalization rate, degradation rate, or
  release yield unless measured or supplied by validated literature.
- MNAzyme sequence efficacy or final cell-killing probability.

## V1.2 restricted extension

V1.2 preserves ten ODE states and compares `linear_v1 + accessibility_v1_1`
for V1.1 compatibility with `saturable_v1_2 + sequential_v1_2` for the V1.2
exploration. It is not a full PBPK model: no lung, kidney, normal-bone,
other-tissue, transit, or vascular states are added. Unmodeled organ passage is
not automatically terminal elimination in this central-circulation approximation.

`C_blood = A_blood / Vblood` is normalized carrier per apparent volume, not an
absolute mass concentration. Fixed `R_EphA2_eff` is not a dynamic receptor or
receptor-saturation model. Model fractions are not probabilities of clinical
success.
