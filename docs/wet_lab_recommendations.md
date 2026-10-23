# Wet-Lab Recommendations

The first executable model should be used to rank measurement priorities, not to
claim validated biological delivery fractions. Current highest-value
measurements are:

1. Carrier size distribution and aggregation state in blood-relevant media.
2. Blood stability or circulation half-life for the actual CaCO3 formulation.
3. Liver and spleen uptake or biodistribution for the actual carrier.
4. CaCO3 degradation rate as a function of pH.
5. Cargo loading and active release yield.
6. OCN-hydroxyapatite retention and release after carrier conjugation.
7. YSA-EphA2 binding and internalization under osteosarcoma-cell conditions.
8. Tumor-model IFP and permeability context for any future quantitative
   biodistribution experiment.

## V1.1 measurement priority

Directly test whether an OCN-bound carrier retains any effective YSA-EphA2
cellular approach, binding, or internalization capacity under relevant matrix
and osteosarcoma-cell conditions. Useful measurements would compare free and
hydroxyapatite/OCN-retained carriers using the same conjugation batch and report
time-resolved binding and internalization. This recommendation identifies an
uncertainty to measure; it does not claim that the modeled OCN-to-YSA transfer
exists or that OCN-YSA must outperform YSA-only.

Interpret state stock, exposure AUC, and cumulative flux separately when
matching experiments to outputs. Paired Monte Carlo fractions are not clinical
success probabilities.

## V1.2 measurement priorities

Measure liver and spleen uptake over a dose range to distinguish a linear
effective uptake approximation from the exploratory saturable-MPS assumption.
Measure whether intact carrier reappears after organ-associated uptake before
attributing any late blood phase to `k_liver_release` or `k_spleen_release`;
MPS release is only one candidate mechanism. In matrix-plus-cell experiments,
test whether OCN-retained double-modified carrier can transition to an
EphA2-associated/internalized population to constrain the sequential
`k_cross` assumption. These experiments should report normalized or absolute
units explicitly; the current `C_blood` is not an absolute mass concentration.

