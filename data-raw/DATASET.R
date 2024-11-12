# CDISC PHUSE PHASE 2 3 armed trial - in ADAM format --------------------

# install.packages("safetyData")
# library(safetyData)
# library(dplyr)

trial02_cdisc <- safetyData::adam_adtte |>
  dplyr::select(
    USUBJID, AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX,
    TRTP, CNSR, SAFFL
  ) |>
  dplyr::rename(
    AVAL = CNSR,
    FASFL = SAFFL
  ) |>
  dplyr::mutate(
    AVAL = 1 - AVAL,
    ## Correct CNSR == 0 is a an event.
    AVALC = dplyr::if_else(AVAL == 1, "YES", "NO"),
    PARAM = "Dermatologic Event (Yes / No)",
    TRTPN = dplyr::case_when(
      TRTP == "Placebo" ~ 1,
      TRTP == "Xanomeline Low Dose" ~ 2,
      TRTP == "Xanomeline High Dose" ~ 3
    )
  ) |>
  dplyr::select(
    USUBJID, PARAM, AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX,
    TRTP, TRTPN, AVAL, AVALC, FASFL
  )


# attach missing var labels
attr(trial02_cdisc$PARAM, "label") <- "Parameter name"
attr(trial02_cdisc$TRTPN, "label") <- "Planned treatment (Numeric)"
attr(trial02_cdisc$AVAL, "label") <- "Analysis value"
attr(trial02_cdisc$AVALC, "label") <- "Analysis value (Character)"
attr(trial02_cdisc$FASFL, "label") <- "Full analysis set flag"


# Hypothetical Phase 3 trial from simulated data ------------------------

# load data - local version
# variable names toLowerCase
trial01 <-
  readr::read_csv("inst/extdata/sim_trial_01.csv",
    show_col_types = FALSE
  ) |>
  dplyr::as_tibble()


# attach var labels
attr(trial01$usubjid, "label") <- "Unique patient identifier"
attr(trial01$aval, "label") <- "Primary outcome variable (yes/no)"
attr(trial01$trtp, "label") <- "Planned treatment"
attr(trial01$bl_cov, "label") <- "Baseline covariate (numeric)"
attr(trial01$bl_cov_c, "label") <- "Baseline covariate (category)"
attr(trial01$region_2, "label") <- "Region 2 indicator (yes/no)"
attr(trial01$region_3, "label") <- "Region 3 indicator (yes/no)"
attr(trial01$region_4, "label") <- "Region 4 indicator (yes/no)"
attr(trial01$region_5, "label") <- "Region 5 indicator (yes/no)"


# SAS %margins macro output for Ge et al. method ------------------------
# on the trial_01 dataset

# load data - local version
margins_trial01 <-
  haven::read_sas("inst/extdata/margins_trt_diffs.sas7bdat")

# SAS ge et al. macro output for Ge et al. method ------------------------
# on the trial_01 dataset

# load data - local version
ge_macro_trial01 <-
  haven::read_sas("inst/extdata/ge_macro_trt_diffs.sas7bdat")

# Load data ------------------------------------------------------
usethis::use_data(trial01, trial02_cdisc, margins_trial01, ge_macro_trial01, overwrite = TRUE, internal = FALSE)
