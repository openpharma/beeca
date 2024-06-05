#' Example CDISC Clinical Trial Dataset in ADaM Format
#'
#' This dataset is a simplified, binary outcome version of a sample Phase 2 clinical trial dataset
#' formatted according to the Analysis Data Model (ADaM) standards set by the Clinical Data
#' Interchange Standards Consortium (CDISC). It is designed for training and educational purposes,
#' showcasing how clinical trial data can be structured for statistical analysis.
#'
#' @format A data frame with 254 rows and 13 columns, representing trial participants and key variables:
#' \describe{
#'   \item{USUBJID}{Unique subject identifier (alphanumeric code). A code unique to the clinical trial}
#'   \item{PARAM}{Parameter name indicating the specific measurement or outcome assessed.}
#'   \item{AGE}{Age of the participant at study enrollment, in years.}
#'   \item{AGEGR1}{Categorical representation of age groups.}
#'   \item{AGEGR1N}{Numeric code representing age groups, used for statistical modeling.}
#'   \item{RACE}{Self-identified race of the participant}
#'   \item{RACEN}{Numeric representation of race categories, used for statistical modeling.}
#'   \item{SEX}{Participant's sex at birth.}
#'   \item{TRTP}{Planned treatment assignment, indicating the specific intervention or control condition.}
#'   \item{TRTPN}{Numeric code for the planned treatment, simplifying data analysis procedures.}
#'   \item{AVAL}{Analysis value, representing the primary outcome measure for each participant.}
#'   \item{AVALC}{Character representation of the analysis value, used in descriptive summaries.}
#'   \item{FASFL}{Full analysis set flag, indicating if the participant's data is included in the full analysis set.}
#' }
#' @source This dataset has been reformatted for educational use from the
#' `safetyData` package, specifically `adam_adtte`. For the original data and
#' more detailed information, please refer to the
#' [`safetyData`](https://safetygraphics.github.io/safetyData/) documentation.
#'
#' @details
#' This dataset serves as an illustrative example for those learning about the
#' ADaM standard in clinical trials. It includes common variables like
#' demographic information, treatment assignments, and outcome measures.
#'
#' Data privacy and ethical considerations have been addressed through the
#' anonymization of subject identifiers and other sensitive information.
#' The dataset is intended for educational and training purposes only.
#'
#' @note
#' The numeric codes for categorical variables such as `RACEN` and `TRTPN`
#' are arbitrary and should be interpreted within the context of this dataset.
#' For example, refer to the categorical representations for additional context.
#'
"trial02_cdisc"
