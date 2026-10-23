# Solver accuracy and historical replay

The first Ubuntu CI run passed 252 tests and failed the strict comparison with
a delivery-model trajectory saved on Windows. Both the failure and its original
threshold (`rtol=1e-9`, `atol=1e-11`) are retained in the numerical record.

Independent comparisons clarified what that check measured:

| Comparison | Maximum absolute difference, normalized units |
| --- | ---: |
| Ubuntu default LSODA vs saved Windows trajectory | 1.219e-9 |
| Default LSODA vs tight Radau | 7.538e-9 |
| Tight LSODA vs tight Radau on Ubuntu | 3.399e-11 |

The saved trajectory was calculated with LSODA at `rtol=1e-8`, `atol=1e-10`.
On Windows, the default result reproduced that trajectory, while more accurate
LSODA and Radau solutions differed from it by about 7.5e-9. Tightening LSODA
tolerances reduced its discrepancy from tight Radau. A many-decimal CSV records
the computed numbers; it does not increase their numerical accuracy.

The checks now answer separate questions:

The seven core model files and parameter CSV also match their frozen SHA256
values. These hashes and the fixture checksum are recorded in
`../research/numerical_reference.json`; numerical comparisons alone would not
establish source identity.

- The ablation test compares its saved time courses with a direct simulation
  on the same numerical backend, retaining the original comparison tolerance.
  This checks that output processing preserves ODE states.
- A solver-convergence test compares default LSODA with tight Radau using a
  normalized whole-trajectory error budget of 1e-7, and compares tight LSODA
  with tight Radau using 1e-9. Mass balance and non-negativity are checked too.
  These are explicit acceptance budgets, not theoretical global-error bounds.
- The original strict frozen-trajectory comparison remains available in the
  diagnostic command. Its Ubuntu result remains a failed historical replay;
  it is not relabeled as a successful replay.

```bash
python scripts/check_solver_agreement.py
python scripts/check_solver_agreement.py --strict-history
```

The first command checks numerical agreement and prints the historical replay
result separately. The second additionally returns a failing exit status if
the original strict historical threshold or recorded source/fixture identity
is not met. It records both checkout-file and actually imported-module hashes,
including the parameter registry actually read. Neither command changes
the model equations, parameters or saved historical fixture.
