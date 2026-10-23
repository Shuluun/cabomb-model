# V1.3 Explicit-Circulation Parameter Evidence

Status: working evidence record for the V1.3 explicit-circulation mode.

This file documents the source basis for every new V1.3 parameter before any
value is committed to `data/processed/parameters_master.csv`. It follows the
applicability order required by `AGENTS.md` and `CLAUDE.md`:

1. CABOMB actual CaCO3 carrier
2. similar CaCO3 carriers
3. similar inorganic nanoparticles
4. generic nanoparticle literature
5. human physiology references (for circulation volumes and flow fractions)

A value that does not come from a CABOMB measurement is marked
`assumed_exploratory` (or `literature_derived` when it rests on a standard
human physiology reference) and is **not** a measured biological value.

## Identifier caveat

Where a journal citation is listed, the bibliographic details are given to
author/year/journal level. Exact DOI/PMID strings are stated only where they
were confirmed in this session's search; where they are not stated, the row
says `confirm against PubMed/CrossRef` and the value is treated as a
transparent reference, not a verified primary record. No citation is invented,
and no docking score is converted into a binding constant.

## Source inventory

| Ref | Source | What it supports |
| --- | --- | --- |
| R1 | Guyton & Hall, *Textbook of Medical Physiology* (standard 70 kg human) | Fractional distribution of blood volume: systemic veins ~64 %, systemic arteries ~13 %, pulmonary ~9 %, heart ~7 %, systemic capillaries ~7 % |
| R2 | Brown R.P., Delp M.D., Lindstedt S.L., Rhomberg L.R., Beliles R.P., *Toxicol. Ind. Health* 1997;13(4):407-84 | Organ blood-flow fractions and tissue blood volumes for PBPK models |
| R3 | Davies B., Morris T., *Pharm. Res.* 1993;10(7):1093-5 | Human organ blood-flow fractions for PBPK models |
| R4 | Kutscher H.L. et al., *J. Control. Release* 2010;143(1):31-7 | Rigid-microparticle pulmonary transit: 2-3 μm pass the lung to liver/spleen; entrapment rises sharply above ~6-7 μm |
| R5 | Project `data/processed/parameters_master.csv` + supplied workbook | CABOMB cardiac output (6.5 L/min), organ flow fractions (liver 22.5 %, spleen 3 %), tumor flow (0.045 L/min), Vblood 5.2 L, Vtum 0.01 L |
| R6 | Project thresholds (`d_lung_mechanical_retention_threshold` = 7 μm, `d_mps_high_uptake_threshold` = 200 nm, `d_tumor_microvessel_pore`) | Size cut-offs used to anchor the size-factor ranges |

The organ **flow fractions already in the registry (R5) are project-supplied**
and are used for V1.3 branch flows; R2/R3 are cited as the physiology
cross-check for the *order* of the fractions, not as a replacement for the
project values. `Q_rest` is derived by mass conservation as
`Qco - (Qliv + Qspl + Qtum)` and is not an independent literature row.

## New V1.3 parameters (34 rows)

Every row below maps to one registry row. `D/E/F` are the registry
`baseline/lower/upper` columns.

### A. Vascular volumes (L)

The nine vascular volumes are allocated so that they are consistent with the
project total `Vblood = 5.2 L` and with the R1 fractional distribution. They
are human-physiology references (not CaCO3-specific). The central-venous
versus rest-vascular venous split is a modelling allocation because the
project supplies only a total blood volume, so those two rows are
`assumed_exploratory`.

| Symbol | D (baseline) | E (lower) | F (upper) | Basis | Classification |
| --- | --- | --- | --- | --- | --- |
| `V_central_venous` | 1.0 | 0.5 | 1.5 | R1 large-vein pool + allocation of venous remainder | assumed_exploratory |
| `V_right_heart` | 0.15 | 0.10 | 0.20 | R1 heart pool split between right/left heart | literature_derived |
| `V_lung_vascular` | 0.45 | 0.35 | 0.55 | R1 pulmonary ~9 % of 5.2 L | literature_derived |
| `V_left_heart` | 0.15 | 0.10 | 0.20 | R1 heart pool split between right/left heart | literature_derived |
| `V_systemic_arterial` | 0.70 | 0.55 | 0.85 | R1 systemic arteries ~13 % of 5.2 L | literature_derived |
| `V_liver_vascular` | 0.10 | 0.05 | 0.20 | R2/R3 hepatic blood volume fraction | assumed_exploratory |
| `V_spleen_vascular` | 0.03 | 0.01 | 0.06 | R2/R3 splenic blood volume fraction | assumed_exploratory |
| `V_tumor_vascular` | 1.5e-4 | 5e-5 | 5e-4 | ~1.5 % of `Vtum = 0.01 L`; mouse tumor vascular fraction ~1-2 % | assumed_exploratory |
| `V_rest_vascular` | 2.62 | 1.8 | 3.4 | R1 residual (small veins + systemic capillaries) | assumed_exploratory |

The nine volumes sum to approximately 5.2 L at baseline. `V_tumor_vascular` is
very small (≈ 0.15 mL), which produces a sub-second tumor-vascular transit and
is the reason the V1.3 solver uses a dense early time grid.

### B. Lung first-pass retention fraction by size (dimensionless)

R4 shows that rigid particles below ~6-7 μm are **not** mechanically entrapped
in the pulmonary capillary bed; 2-3 μm particles transit the lung and localize
in liver/spleen. The project threshold (R6) sets mechanical retention at
>7 μm, outside every CABOMB size scenario. The fractions below therefore
represent only the small pulmonary-macrophage / transient-sequestration
component, not mechanical filtration. All are `assumed_exploratory`, low
confidence.

| Symbol | D | E | F |
| --- | --- | --- | --- |
| `f_lung_ret_100nm` | 0.002 | 0.0 | 0.05 |
| `f_lung_ret_200nm` | 0.002 | 0.0 | 0.05 |
| `f_lung_ret_500nm` | 0.005 | 0.0 | 0.08 |
| `f_lung_ret_800nm` | 0.01 | 0.0 | 0.10 |
| `f_lung_ret_2um` | 0.015 | 0.0 | 0.05 |
| `f_lung_ret_3um` | 0.025 | 0.0 | 0.07 |

### C. MPS (liver/spleen) size factor (dimensionless)

R6 `d_mps_high_uptake_threshold = 200 nm` states that carriers above ~200 nm
are increasingly taken up by the RES/MPS. R4 shows 2-3 μm particles localize
80-90 % to liver/spleen. The factor multiplies `Vmax_liver`/`Vmax_spleen`.
All `assumed_exploratory`, low confidence.

| Symbol | D | E | F |
| --- | --- | --- | --- |
| `size_mps_100nm` | 1.0 | 0.8 | 1.5 |
| `size_mps_200nm` | 1.3 | 0.9 | 2.0 |
| `size_mps_500nm` | 1.8 | 1.2 | 3.0 |
| `size_mps_800nm` | 2.2 | 1.4 | 4.0 |
| `size_mps_2um` | 2.8 | 1.7 | 4.5 |
| `size_mps_3um` | 3.2 | 2.0 | 5.5 |

### D. Tumor extravasation size factor (dimensionless)

R6 `d_tumor_microvessel_pore` sets the tumor vessel pore window; sub-200 nm
carriers extravasate well, 500-800 nm poorly, and 2-3 μm are essentially
excluded from the EPR pore. The factor multiplies `k_entry_eff`. All
`assumed_exploratory`, low confidence.

| Symbol | D | E | F |
| --- | --- | --- | --- |
| `size_extra_100nm` | 1.0 | 0.8 | 1.2 |
| `size_extra_200nm` | 0.9 | 0.7 | 1.1 |
| `size_extra_500nm` | 0.4 | 0.2 | 0.6 |
| `size_extra_800nm` | 0.2 | 0.05 | 0.4 |
| `size_extra_2um` | 0.03 | 0.0 | 0.12 |
| `size_extra_3um` | 0.015 | 0.0 | 0.08 |

### E. Internalization size factor (dimensionless)

Small (100-200 nm) carriers internalize most efficiently; larger carriers
internalize progressively less. The factor multiplies both `k_int` and
`k_int_ns`. All `assumed_exploratory`, low confidence.

| Symbol | D | E | F |
| --- | --- | --- | --- |
| `size_int_100nm` | 1.0 | 0.8 | 1.2 |
| `size_int_200nm` | 1.0 | 0.8 | 1.2 |
| `size_int_500nm` | 0.6 | 0.3 | 1.0 |
| `size_int_800nm` | 0.4 | 0.2 | 0.8 |
| `size_int_2um` | 0.12 | 0.0 | 0.32 |
| `size_int_3um` | 0.08 | 0.0 | 0.25 |

### F. Direct free-to-YSA accessibility (dimensionless)

`alpha_ysa_direct_access` gates the direct free-to-YSA route in
`parallel_dual_access_v1_3` only. It is a structural accessibility assumption
for the parallel dual-access hypothesis, not a measured binding constant.

| Symbol | D | E | F |
| --- | --- | --- | --- |
| `alpha_ysa_direct_access` | 0.5 | 0.0 | 1.0 |

`0` recovers the sequential route exactly; `1` gives full direct free-to-YSA
access. All `assumed_exploratory`, low confidence.

## Size scenario labels

| Label | `size` key | Meaning |
| --- | --- | --- |
| `100nm` | `100nm` | original design scenario |
| `200nm` | `200nm` | optimization target |
| `500nm` | `500nm` | oral/EM observation, full-distribution exploratory |
| `800nm` | `800nm` | oral/EM observation, full-distribution exploratory |
| `2um` | `2um` | oral/EM observation, full-distribution exploratory |
| `3um` | `3um` | oral/EM observation, full-distribution exploratory |

100 nm and 200 nm are full-distribution scenarios; 500/800 nm, 2 μm, and 3 μm
are treated as exploratory extensions per the project brief, not as validated
distribution data.

## Explicit gaps

- No CABOMB-measured vascular volume, transit time, or size-factor value was
  supplied. All V1.3 rows are transparent physiology references or exploratory
  assumptions and must not be read as measured.
- The central-venous / rest-vascular venous allocation is a structural
  modelling choice, not a measurement.
- Size factors are monotone, order-preserving extrapolations anchored to the
  project size thresholds (R6); they are not calibrated against any CaCO3
  biodistribution dataset.
