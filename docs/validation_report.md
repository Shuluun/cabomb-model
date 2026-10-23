# Validation of the release copy

Checked locally on 2026-10-23, Windows, Python 3.12.14.

| Check | Result |
| --- | --- |
| Wheel build and installation into a fresh environment | Passed |
| Installed source and resources match this release copy | Passed |
| Complete test suite after numerical portability repair | 254 passed, 462.36 seconds |
| Parameter registry validation | Passed |
| V1.3 baseline and full report generation | Passed |
| Synthetic early-only forecast and separate late scoring | Passed |
| Fixed-input intracellular response analysis | All numerical checks passed |
| Conditional calcium-equivalent calculator | Passed |
| Streamlit app execution, default and V1.3 selection | No execution exceptions |
| Local headless Streamlit service | HTTP 200 health response |
| Installed dependency consistency | `pip check` passed |
| Reviewed analysis scope | Five active grid dimensions; inactive sequential routes excluded |

The earlier clean-wheel check passed 253 tests. The latest release-copy recheck
includes the added independent-solver convergence test.

The later Ubuntu solver-history comparison and its numerical diagnosis are
documented in [Numerical portability](numerical_portability.md). Output processing,
solver convergence and strict historical replay are evaluated separately.

The model and forecast commands ran from a directory outside the public source
tree using the wheel-installed modules. Parameter resources were loaded from
that installation. The original project was not used as a runtime fallback.

The bridge's largest first-order reference error was 7.42666537e-08 normalized
units; the largest mid/fine quadrature difference was 1.14591861e-07. Synthetic
and bridge figures were inspected for labels and layout.

Tested versions: cabomb-model 0.2.0, numpy 2.5.4, scipy 1.18.1, pandas 3.1.0, matplotlib 3.11.2, openpyxl 3.1.5, pytest 9.1.1, streamlit 1.66.1.
`requirements-tested.txt` records the environment. To reproduce it, install that
file and then install the project. The main README uses supported dependency
ranges for a simpler setup.

The table above records the clean Windows installation. GitHub workflow runs
are listed in the repository's [Actions tab](https://github.com/Shuluun/cabomb-model/actions). The Streamlit version
emits a deprecation notice for an existing width option, without affecting the
tested default and V1.3 views.

These checks establish software execution and numerical behavior. The prior
literature comparison and its separate historical replay check are described
in `research/findings.md`.
