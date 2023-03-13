# Name: ADAE
#
# Label: Adverse Event Analysis Dataset
#
# Input: ae, adsl, ex_single
library(admiral)
library(dplyr)
library(lubridate)
library(stringr)

# Load source datasets ----

ae <- convert_blanks_to_na(read_xpt("sdtm/ae.xpt"))
suppae <- convert_blanks_to_na(read_xpt("sdtm/suppae.xpt"))

adsl <- convert_blanks_to_na(read_xpt("adam/adsl.xpt"))

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(USUBJID, STUDYID, RACE, RACEN, SAFFL, SEX, SITEID, TRT01A,
                  TRT01AN, TRTSDT, TRTEDT, AGE, AGEGR1, AGEGR1N)

adae0 <- ae %>%
  # join adsl to ae
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by = vars(STUDYID, USUBJID)
  ) %>%
  rename(TRTA=TRT01A, TRTAN=TRT01AN) %>%
  ## Derive analysis start time ----
  derive_vars_dt(
    dtc = AESTDTC,
    new_vars_prefix = "AST",
    highest_imputation = "D",
    date_imputation = "first",
    flag_imputation = "auto",
    min_dates = vars(TRTSDT)
  ) %>%
  ## Derive analysis end time ----
  derive_vars_dt(
    dtc = AEENDTC,
    new_vars_prefix = "AEN",
    highest_imputation = "n"
  ) %>%
  ## Derive analysis duration (value and unit) ----
  restrict_derivation(
    derivation = derive_vars_duration,
    args = params(
      new_var = ADURN,
      new_var_unit = ADURU,
      start_date = ASTDT,
      end_date = AENDT,
      in_unit = "days",
      out_unit = "days",
      add_one = TRUE,
      trunc_out = FALSE
    ),
    filter = is.na(ASTDTF)
  ) %>%
  mutate(ADURU=case_when(ADURU == 'DAYS' ~ 'DAY', T ~ ADURU)) %>%
  ## Derive analysis start relative day and  analysis end relative day ----
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = vars(ASTDT, AENDT)
  ) %>%
  ## Derive CQ01NAM ----
  mutate(CQ01NAM=case_when(
    str_detect(str_to_upper(AEDECOD), 'APPLICATION') ~ 'DERMATOLOGIC EVENTS',
    str_detect(str_to_upper(AEDECOD), 'DERMATITIS') ~ 'DERMATOLOGIC EVENTS',
    str_detect(str_to_upper(AEDECOD), 'ERYTHEMA') ~ 'DERMATOLOGIC EVENTS',
    str_detect(str_to_upper(AEDECOD), 'BLISTER') ~ 'DERMATOLOGIC EVENTS',
    AEBODSYS == 'SKIN AND SUBCUTANEOUS TISSUE DISORDERS'
      & !(str_to_upper(AEDECOD) %in% c('COLD SWEAT', 'HYPERHIDROSIS', 'ALOPECIA'))
      ~ 'DERMATOLOGIC EVENTS',
    TRUE ~ NA_character_
  )) %>%
  ## Derive TRTEMFL ----
  mutate(TRTEMFL=case_when(ASTDT >= TRTSDT ~ 'Y', T ~ NA_character_))

  ## Derive occurrence flags ----
adae1 <- adae0 %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID),
      order = vars(USUBJID, ASTDT, AESEQ),
      new_var = AOCC01FL,
      mode = "first"
    ),
    filter = CQ01NAM == 'DERMATOLOGIC EVENTS' & TRTEMFL == 'Y'
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID),
      order = vars(USUBJID, ASTDT, AESEQ),
      new_var = AOCC02FL,
      mode = "first"
    ),
    filter = AESER == 'Y' & TRTEMFL == 'Y'
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID, AEBODSYS),
      order = vars(USUBJID, AEBODSYS, ASTDT, AESEQ),
      new_var = AOCC03FL,
      mode = "first"
    ),
    filter = AESER == 'Y' & TRTEMFL == 'Y'
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID, AEBODSYS, AEDECOD),
      order = vars(USUBJID, AEBODSYS, AEDECOD, ASTDT, AESEQ),
      new_var = AOCC04FL,
      mode = "first"
    ),
    filter = AESER == 'Y' & TRTEMFL == 'Y'
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID),
      order = vars(USUBJID, ASTDT, AESEQ),
      new_var = AOCCFL,
      mode = "first"
    ),
    filter = TRTEMFL == 'Y'
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID, AEBODSYS, AEDECOD),
      order = vars(USUBJID, AEBODSYS, AEDECOD, ASTDT, AESEQ),
      new_var = AOCCPFL,
      mode = "first"
    ),
    filter = TRTEMFL == 'Y'
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID, AEBODSYS),
      order = vars(USUBJID, AEBODSYS, ASTDT, AESEQ),
      new_var = AOCCSFL,
      mode = "first"
    ),
    filter = TRTEMFL == 'Y'
  ) %>%
  # Final select ----
  select(STUDYID, SITEID, USUBJID, TRTA, TRTAN, AGE, AGEGR1, AGEGR1N, RACE,
         RACEN, SEX, SAFFL, TRTSDT, TRTEDT, ASTDT, ASTDTF, ASTDY, AENDT, AENDY,
         ADURN, ADURU, AETERM, AELLT, AELLTCD, AEDECOD, AEPTCD, AEHLT, AEHLTCD,
         AEHLGT, AEHLGTCD, AEBODSYS, AESOC, AESOCCD, AESEV, AESER, AESCAN,
         AESCONG, AESDISAB, AESDTH, AESHOSP, AESLIFE, AESOD, AEREL, AEACN,
         AEOUT, AESEQ, TRTEMFL, AOCCFL, AOCCSFL, AOCCPFL, AOCC02FL, AOCC03FL,
         AOCC04FL, CQ01NAM, AOCC01FL)

# Add labels

adae_labels <- c(
  'Study Identifier', 'Study Site Identifier', 'Unique Subject Identifier',
  'Actual Treatment', 'Actual Treatment (N)', 'Age', 'Pooled Age Group 1',
  'Pooled Age Group 1 (N)', 'Race', 'Race (N)', 'Sex', 'Safety Population Flag',
  'Date of First Exposure to Treatment', 'Date of Last Exposure to Treatment',
  'Analysis Start Date', 'Analysis Start Date Imputation Flag',
  'Analysis Start Relative Day', 'Analysis End Date',
  'Analysis End Relative Day', 'AE Duration (N)', 'AE Duration Units',
  'Reported Term for the Adverse Event', 'Lowest Level Term',
  'Lowest Level Term Code', 'Dictionary-Derived Term', 'Preferred Term Code',
  'High Level Term', 'High Level Term Code', 'High Level Group Term',
  'High Level Group Term Code', 'Body System or Organ Class',
  'Primary System Organ Class', 'Primary System Organ Class Code',
  'Severity/Intensity', 'Serious Event', 'Involves Cancer',
  'Congenital Anomaly or Birth Defect',
  'Persist or Signif Disability/Incapacity', 'Results in Death',
  'Requires or Prolongs Hospitalization', 'Is Life Threatening',
  'Occurred with Overdose', 'Causality', 'Action Taken with Study Treatment',
  'Outcome of Adverse Event', 'Sequence Number',
  'Treatment Emergent Analysis Flag', '1st Occurrence of Any AE Flag',
  '1st Occurrence of SOC Flag', '1st Occurrence of Preferred Term Flag',
  '1st Occurrence 02 Flag for Serious',
  '1st Occurrence 03 Flag for Serious SOC',
  '1st Occurrence 04 Flag for Serious PT', 'Customized Query 01 Name',
  '1st Occurrence 01 Flag for CQ01'
)

adae2 <- data.frame(mapply(
  function (var, lab) {
    col <- adae1[[var]]
    attr(col, "label") <- lab
    return(col)
  },
  var=colnames(adae1), lab=adae_labels,
  SIMPLIFY = FALSE
))

# SAS format attributes

attr(adae2$AENDT, "format.sas") <- "DATE9"
attr(adae2$ASTDT, "format.sas") <- "DATE9"
attr(adae2$TRTEDT, "format.sas") <- "DATE9"
attr(adae2$TRTSDT, "format.sas") <- "DATE9"

# Save output ----

xportr_write(adae2, "adam/adae.xpt")
