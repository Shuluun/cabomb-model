# Open Questions

These questions are prioritized by expected impact on first-version model
conclusions and by whether the missing information can change the design ranking.

## High Priority

1. What is the actual CABOMB CaCO3 carrier size distribution in blood-relevant
   media?
   - Impact: particle size affects MPS uptake, lung filtration risk, tumor
     extravasation, interstitial diffusion, and internalization.
   - Current status: project files mention 100 nm context and broad exploratory
     particle-size scans, but no final measured distribution is supplied.

2. What are the CaCO3 carrier blood stability and circulation half-life?
   - Impact: this strongly controls how much carrier remains available for tumor
     entry.
   - Current status: project files identify half-life as missing and recommend
     scans, not a measured value.

3. What are liver and spleen uptake/clearance rates for the actual CaCO3
   formulation?
   - Impact: MPS uptake may dominate biodistribution.
   - Current status: spreadsheet gives broad project reference values for liver
     uptake, but no CaCO3-specific experiment or citation.

4. What is tumor entry/permeability for osteosarcoma-relevant carriers?
   - Impact: low vs high EPR may reverse whether targeting differences matter.
   - Current status: literature says EPR is heterogeneous; no osteosarcoma
     CaCO3-specific permeability measurement is supplied.

5. What is the actual IFP range in the intended osteosarcoma setting?
   - Impact: IFP directly opposes transvascular/interstitial transport.
   - Current status: spreadsheet gives a human tumor range, but no tumor-type or
     source citation.

6. What are YSA-EphA2 binding parameters under relevant conditions?
   - Impact: needed for YSA-only and OCN-YSA comparisons.
   - Current status: project files ask another teammate for binding-energy data.
     Docking scores must not be treated as experimental `Kd`.

7. What is OCN affinity/retention on hydroxyapatite or bone matrix after
   conjugation to the carrier?
   - Impact: needed to quantify OCN-only and OCN-YSA retention.
   - Current status: project PPT states the qualitative interaction, but no
     quantitative binding or retention data are supplied.

8. What is CaCO3 degradation rate as a function of pH for the actual carrier?
   - Impact: needed to compare premature blood degradation vs intracellular
     release.
   - Current status: project files state acid-dependent degradation but provide
     no rate curve.

## Medium Priority

1. What is EphA2 expression density in selected osteosarcoma cells and relevant
   normal cells?
   - Impact: needed to convert YSA binding into selectivity and internalization
     estimates.

2. What is the receptor-mediated internalization rate after YSA-EphA2 binding?
   - Impact: affects intracellular carrier and cargo-release timing.

3. What dose unit should the model use for the first executable version:
   particle count, CaCO3 mass, cargo amount, or a paired carrier/cargo unit?
   - Impact: affects mass balance, cargo release, and wet-lab interpretability.

4. What is the tumor compartment volume or reference tumor size for simulations?
   - Impact: required for concentration outputs and concentration-dependent
     binding if used.

5. What are cargo loading and release yield per carrier?
   - Impact: needed to connect carrier degradation to released MNAzyme, Fe3+,
     and Ca2+ amounts.

6. Is the intended route strictly IV bolus, infusion, or another injection
   protocol?
   - Impact: initial conditions and input functions differ.

7. Which species should be the first quantitative reference: human exploratory,
   mouse xenograft, or in-vitro-only wet-lab design?
   - Impact: organ volumes, blood flow, EPR, MPS, and scaling all depend on
     species.

8. How should the model treat the project conflict about active targeting and
   `F_target`/`F_extra`?
   - Impact: PPT slide 6 embedded image says active targeting may increase these
     quantities, while slide text and AGENTS.md say it should not.
   - Proposed resolution: first version follows AGENTS.md and slide text:
     active targeting affects binding, retention, and internalization only.

## Low Priority

1. Should lung be an explicit compartment in version 1?
   - Current recommendation: not initially, unless carrier size/aggregation
     suggests lung filtration is material.

2. Should renal clearance be represented for intact carriers?
   - Current recommendation: likely not for 100 nm-scale carriers, but include
     only if measurements indicate small degradation products or free cargo need
     renal handling.

3. Should ECM/collagen density have an explicit parameter?
   - Current recommendation: fold it into low/medium/high EPR or tumor-entry
     scenarios for version 1.

4. Should receptor recycling be modeled?
   - Current recommendation: omit until YSA-EphA2 internalization data exist.

5. Should MNAzyme sequence design be coupled to the delivery model?
   - Current recommendation: keep separate; delivery model can output released
     cargo timing, but not sequence efficacy.

6. Should wet-lab recommendations use absolute concentrations?
   - Current recommendation: only after dose, volume, half-life, uptake,
     degradation, and cargo-yield parameters are sourced or measured.

## Freeze Review Classification

### Questions that can change the model definition

These must be decided before treating the specification as fully frozen for
implementation:

1. Dose and state unit choice: particle count, CaCO3 mass, cargo amount, or a
   paired carrier/cargo unit.
   - Reason: this determines state units, conservation tests, cargo yield
     interpretation, and wet-lab output units.

2. First quantitative reference system: human exploratory, mouse xenograft, or
   in-vitro-only wet-lab design.
   - Reason: compartment volumes, blood flow, MPS clearance, tumor entry, and
     scaling assumptions are species- and system-dependent.

3. Actual carrier size distribution and aggregation state in blood-relevant
   media.
   - Reason: a large or aggregated carrier could require an explicit lung
     filtration compartment or change the validity of the minimum PBPK-inspired
     structure.

4. Administration protocol: IV bolus versus infusion.
   - Reason: bolus uses initial condition `A_blood(0) = D0`; infusion requires
     an input function.

5. Resolution of active-targeting conflict around `F_target`/`F_extra`.
   - Reason: if the project chooses to let active targeting alter tumor entry,
     the current model definition would be invalid. The current freeze follows
     AGENTS.md: active targeting affects only local binding, retention, and
     internalization after tumor-region arrival.

6. Whether cargo accounting must represent only active intracellular release or
   also extracellular/systemic release.
   - Reason: the current specification treats non-intracellular terminal carrier
     loss as lost or inactive payload. Evidence for active extracellular release
     would require additional cargo states.

### Missing parameters that can proceed with explicit exploratory ranges

These should be marked as missing or assumed in the parameter registry, assigned
transparent exploratory ranges, and subjected to sensitivity/uncertainty
analysis before making qualitative design comparisons:

1. Blood circulation half-life and `k_blood_clear`.
2. Liver and spleen uptake/clearance rates for the actual CABOMB formulation.
3. Tumor entry function parameters, including EPR scenario factor,
   permeability, IFP effect, and tumor blood-flow context.
4. OCN-hydroxyapatite association, dissociation, and OCN-bound loss rates.
5. YSA-EphA2 association, dissociation, receptor availability, and
   receptor-mediated internalization rates.
6. pH-dependent CaCO3 degradation rates at physiological and acidic pH.
7. Tumor compartment volume or reference tumor size, if only amount-based
   outputs are interpreted cautiously until concentration conversion is needed.
8. Cargo loading/release yield and active cargo loss rate, if outputs are
   reported as exploratory released-cargo indices rather than quantitative
   therapeutic concentrations.

## V1.2 open questions

1. Does CABOMB carrier uptake by liver and spleen saturate over the intended
   dose range, and what apparent `Vmax`/`Km` values are supported by data?
2. Is organ-associated late loss intact carrier release, degradation, or another
   unmodeled process? MPS release is only one candidate mechanism.
3. Does the double-modified carrier follow a strict OCN-first, YSA-second route,
   or is the sequential structure assumption inadequate?
4. What receptor measurements are required before replacing fixed
   `R_EphA2_eff` with a dynamic receptor/saturation model?
5. Which unmodeled-organ data would justify expanding the central circulation
   approximation rather than treating passage as nonterminal by default?
