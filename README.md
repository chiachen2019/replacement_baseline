# Baseline mosquito abundance and Wolbachia replacement: simulation code

R code accompanying the manuscript **“Baseline mosquito abundance determines minimum viral blocking effect required for Wolbachia replacement as a dengue control tool”** (Chang et al.November 2026).

## Overview

This repository contains deterministic mosquito population simulations investigating how initial *Aedes* abundance, existing source reduction, and Wolbachia-associated viral blocking influence **vectorial capacity**. The model represents aquatic and adult wild-type and Wolbachia-infected mosquitoes, density dependence, seasonality, Wolbachia fitness costs, and discrete releases. Vectorial capacity is calculated from the simulated adult populations.

Two main settings are evaluated:

1. **Existing source reduction:** Vary pre-release reductions in aquatic carrying capacity, the fraction of source reduction retained after Wolbachia introduction, and viral blocking efficacy.
2. **Initially low mosquito abundance:** Vary the initial wild-type population relative to equilibrium and compare Wolbachia replacement with no release and source reduction alone.

The repository contains the simulation and vectorial-capacity calculation scripts, **not** the manuscript's figure-generation scripts or the supplementary analytical derivation.

## Repository structure

```text
.
├── README.md
├── run_simulation.R                   # Main script: runs all four groups of scenarios
├── R/
│   ├── parameters/
│   │   ├── vector_control.R           # Initial conditions and parameters with source reduction
│   │   └── no_vector_control.R        # Initial conditions and parameters without source reduction
│   ├── vector_control_model.R         # Four-compartment ODE model
│   ├── calculate_vectorial_capacity.R # Vectorial-capacity calculation and normalization
│   ├── run_k0_blocking_scenarios.R    # Existing-source-reduction scenarios
│   └── run_initial_blocking_scenarios.R # Initial-population scenarios
└── output/                            # Generated CSV files 
```

## Requirements

- R with native pipe support (`|>`, R 4.1 or later)
- R packages: `deSolve`, `ggplot2`, `tidyverse`

Install packages in R:

```r
install.packages(c("deSolve", "ggplot2", "tidyverse"))
```

The ODE system is solved using `deSolve::ode()`. The main script also loads `ggplot2` and `tidyverse`; `dplyr::bind_rows()` is used to combine results.

## Running the simulations

Clone or download the repository, open a terminal **in the repository root**, and run:

```bash
Rscript run_simulation.R
```

Alternatively, open `run_simulation.R` in RStudio, set the working directory to the repository root, and source the script. The run can be computationally intensive because it integrates many scenarios over a 15-year daily time grid.

The script generates four date-stamped CSV files under `output/`:

| Output prefix | Contents |
| --- | --- |
| `output_partial_` | Wolbachia replacement with varying baseline source reduction, post-release retention, and blocking efficacy |
| `output_initial_blocking_results_` | Wolbachia replacement with varying initial population proportions and blocking efficacy |
| `no_release_results_` | No-Wolbachia-release comparison across initial population proportions |
| `vector_control_only_results_` | Source reduction without Wolbachia releases, across reduction intensities and initial population proportions |

Output files contain simulated time series and scenario identifiers. Relevant columns include `time`, `A_u`, `A_w`, `VC_u`, `VC_w`, `VC_total`, and `VC_relative` (where requested). `VC_relative` is the ratio to the selected baseline, so values below 1 indicate reduced vectorial capacity.

**Baseline definitions:** The existing-source-reduction analysis uses the mean vectorial capacity during the pre-release source-reduction period (`baseline_type = "control_period"`). The initial-population analysis uses vectorial capacity at the initial simulation time (`baseline_type = "initial"`).

## Model and scenarios

The state variables are `L_u` (wild-type aquatic stage), `A_u` (wild-type adults), `L_w` (Wolbachia-infected aquatic stage), and `A_w` (Wolbachia-infected adults). The baseline carrying capacity is `K0 = 30000`, with seasonal amplitude `eta = 0.2` and period 365 days. Releases occur at 7-day intervals for one year; release sizes are calculated by the scenario functions relative to adult wild-type abundance.

The main script evaluates:

- **Existing source reduction:** 9 levels (0.1–0.9) × 4 post-release retention levels (0, 0.25, 0.5, 0.75) × 6 blocking levels (0.5–1.0): **216 combinations**.
- **Initially low abundance, Wolbachia replacement:** 10 initial population fractions (0.1–1.0) × 6 blocking levels: **60 combinations**.
- **No release:** 10 initial population fractions, blocking set to zero, release ratio set to zero.
- **Source reduction only:** 9 source-reduction levels × 11 initial population fractions (0.1–1.0), no Wolbachia release.

See the parameter files and `run_simulation.R` for exact settings. Model parameters and vectorial-capacity defaults are defined in `R/parameters/` and `R/calculate_vectorial_capacity.R`.

## Reproducibility and scope

The source files were reorganized and renamed for portability. The absolute Windows working-directory path was removed, `source()` paths were updated, and generated CSV files were directed to `output/`. The scientific functions and scenario parameter values were otherwise retained as supplied. No numerical equivalence check or full execution has yet been performed for this reorganized version.

**Before citing the repository as an exact reproduction of published results**, verify the code against the final manuscript and archived outputs. In particular, the supplied manuscript describes **54** initial-population Wolbachia scenarios whereas the supplied main script specifies **60**; and the manuscript's reported initial equilibrium populations differ slightly from the hard-coded initial state in `vector_control.R`. These differences have not been silently changed here.

## Citation

Please cite the accompanying manuscript (add journal, year, DOI and final author list when available). For software citation, cite this repository's permanent release/DOI once assigned.

## License

No license is included. Repository maintainers should select and add a license (for example, MIT) after confirming authorization from all code contributors and relevant institutions.

## Contact

Please use the repository's GitHub Issues page for questions or reproducibility reports once the repository is published.
