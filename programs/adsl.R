library(haven)
library(admiral)
library(dplyr)
library(tidyr)
library(metacore)
library(metatools)
library(xportr)

# Data ----

dm <- convert_blanks_to_na(read_xpt("sdtm/dm.xpt"))
ds <- convert_blanks_to_na(read_xpt("sdtm/ds.xpt"))
ex <- convert_blanks_to_na(read_xpt("sdtm/ex.xpt"))
mh <- convert_blanks_to_na(read_xpt("sdtm/mh.xpt"))
sc <- convert_blanks_to_na(read_xpt("sdtm/sc.xpt"))
sv <- convert_blanks_to_na(read_xpt("sdtm/sv.xpt"))
qs <- convert_blanks_to_na(read_xpt("sdtm/qs.xpt"))
ts <- convert_blanks_to_na(read_xpt("sdtm/ts.xpt"))
vs <- convert_blanks_to_na(read_xpt("sdtm/vs.xpt"))

# Grouping functions ---

format_agegr1n <- function(x) {
  case_when(
    x < 65 ~ 1,
    between(x, 65, 80) ~ 2,
    x > 80 ~ 3,
    TRUE ~ NA_real_
  )
}

format_agegr1 <- function(x) {
  case_when(
    x == 1 ~ "<65",
    x == 2 ~ "65-80",
    x == 3 ~ ">80",
    TRUE ~ NA_character_
  )
}

format_disccd <- function(x) {
  case_when(
    x == "SCREEN FAILURE" ~ NA_character_,
    TRUE ~ x
  )
}

format_discreas <- function(x, y) {
  case_when(
    x == "PROTOCOL ENTRY CRITERIA NOT MET" & y == "PROTOCOL VIOLATION" ~ "I/E Not Met",
    x == "PROTOCOL ENTRY CRITERIA NOT MET" ~ NA_character_,
    y == "ADVERSE EVENT" ~ "Adverse Event",
    y == "STUDY TERMINATED BY SPONSOR" ~ "Sponsor Decision",
    y == "DEATH" ~ "Death",
    y == "WITHDRAWAL BY SUBJECT" ~ "Withdrew Consent",
    y == "PHYSICIAN DECISION" ~ "Physician Decision",
    y == "PROTOCOL VIOLATION" ~ "Protocol Violation",
    y == "LOST TO FOLLOW-UP" ~ "Lost to Follow-up",
    y == "LACK OF EFFICACY" ~ "Lack of Efficacy",
    TRUE ~ NA_character_
  )
}

format_arm <- function(x) {
  case_when(
    x == "Screen Failure" ~ NA_character_,
    TRUE ~ x
  )
}

format_armn  <- function(x) {
  case_when(
    x == "Placebo" ~ 0,
    x == "Xanomeline Low Dose" ~ 54,
    x == "Xanomeline High Dose" ~ 81,
    TRUE ~ NA_real_
  )
}

format_eosstt <- function(x) {
  case_when(
    x == "COMPLETED" ~ "COMPLETED",
    TRUE ~ "DISCONTINUED"
  )
}

format_racen <- function(x) {
  case_when(
    x == "AMERICAN INDIAN OR ALASKA NATIVE" ~ 1,
    x == "ASIAN" ~ 2,
    x == "BLACK OR AFRICAN AMERICAN" ~ 3,
    x == "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" ~ 5,
    x == "WHITE" ~ 6,
    TRUE ~ NA_real_
  )
}

format_bmicat <- function(x) {
  case_when(
    x < 25 ~ "<25",
    x < 30 ~ "25-<30",
    x >= 30 ~ ">=30",
    TRUE ~ NA_character_
  )
}

format_durdisc <- function(x) {
  case_when(
    x < 12 ~ "<12",
    x >= 12 ~ ">=12",
    TRUE ~ NA_character_
  )
}


format_sitegr1 <- function(x) {
  return (case_when(
    x %in% c("715", "717") ~ "900",
    TRUE ~ x
  ))
}

# Derivations on source domains ----

adsl_dm <- dm %>%
  select(STUDYID, USUBJID, SUBJID, SITEID, SEX, RFSTDTC, RFENDTC, RACE, ETHNIC,
         AGE, AGEU, ARM, ARMCD, DTHFL, ACTARM) %>%
  mutate(AGEGR1N = format_agegr1n(AGE)) %>%
  mutate(AGEGR1 = format_agegr1(AGEGR1N)) %>%
  mutate(RACEN = format_racen(RACE)) %>%
  mutate(ARMN = format_armn(ARM)) %>%
  mutate(TRT01P = ARM) %>%
  mutate(TRT01A = ARM) %>%
  mutate(TRT01PN = format_armn(TRT01P)) %>%
  mutate(TRT01AN = format_armn(TRT01A)) %>%
  mutate(SITEGR1 = format_sitegr1(SITEID)) %>%
  mutate(ITTFL = case_when(ARMCD == "Scrnfail" ~ "N", T ~ "Y")) %>%
  derive_vars_dt(
    dtc = RFENDTC,
    new_vars_prefix = "RFEN",
    flag_imputation = "none"
  )

adsl_sv <- sv %>%
  select(USUBJID, VISITNUM, SVSTDTC) %>%
  filter(VISITNUM %in% c(1, 3, 4, 8, 10, 12)) %>%
  pivot_wider(id_cols=USUBJID, names_from=VISITNUM, values_from=SVSTDTC, names_prefix="VISIT") %>%
  derive_vars_dt(
    dtc = VISIT3,
    new_vars_prefix = "TRTS",
    flag_imputation = "none"
  ) %>%
  derive_vars_dt(
    dtc = VISIT1,
    new_vars_prefix = "VISIT1",
    flag_imputation = "none"
  ) %>%
  derive_vars_dt(
    dtc = VISIT4,
    new_vars_prefix = "VISIT4",
    flag_imputation = "none"
  ) %>%
  derive_vars_dt(
    dtc = VISIT8,
    new_vars_prefix = "VISIT8",
    flag_imputation = "none"
  ) %>%
  derive_vars_dt(
    dtc = VISIT10,
    new_vars_prefix = "VISIT10",
    flag_imputation = "none"
  ) %>%
  derive_vars_dt(
    dtc = VISIT12,
    new_vars_prefix = "VISIT12",
    flag_imputation = "none"
  ) %>%
  select(-c(VISIT1, VISIT3, VISIT4, VISIT8, VISIT10, VISIT12))

adsl_ex <- ex %>%
  select(USUBJID, EXENDTC, EXSEQ) %>%
  arrange(USUBJID, EXSEQ) %>%
  group_by(USUBJID) %>%
  filter(EXSEQ == last(EXSEQ)) %>%
  ungroup() %>%
  select(-EXSEQ) %>%
  derive_vars_dt(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    flag_imputation = "none"
  ) %>%
  select(-EXENDTC)

adsl_ds_disp <- ds %>%
  filter(DSCAT == "DISPOSITION EVENT") %>%
  select(USUBJID, DSTERM, DSDECOD, VISITNUM, DSDTC) %>%
  mutate(DCDECOD = format_disccd(DSDECOD)) %>%
  mutate(DCSREAS = format_discreas(DSTERM, DCDECOD)) %>%
  mutate(DISCONFL = case_when(DCSREAS != "COMPLETED" & !is.na(DCSREAS) ~ "Y", T ~ NA_character_)) %>%
  mutate(DSRAEFL = case_when((DCSREAS) == "Adverse Event" ~ "Y", T ~ NA_character_)) %>%
  derive_vars_dt(
    dtc = DSDTC,
    new_vars_prefix = "DS",
    flag_imputation = "none"
  ) %>%
  select(-DSDTC) %>%
  rename(DSVSTNUM=VISITNUM)

adsl_qs <- qs %>%
  filter(QSCAT == "MINI-MENTAL STATE") %>%
  select(USUBJID, QSORRES) %>%
  group_by(USUBJID) %>%
  summarize(MMSETOT=sum(as.numeric(QSORRES)), .groups="drop")

adsl0 <- adsl_dm %>%
  filter(ACTARM != "Screen Failure") %>%
  left_join(adsl_ds_disp, by = "USUBJID") %>%
  left_join(adsl_ex, by = "USUBJID") %>%
  left_join(adsl_qs, by = "USUBJID") %>%
  left_join(adsl_sv, by = "USUBJID")

# ADSL derivations ----

adsl1 <- adsl0 %>%
  ## Derive TRTEDT ----
  mutate(TRTEDT=case_when(
    is.na(EXENDT) & DSVSTNUM > 3 & DSDECOD != "COMPLETED" ~ DSDT,
    TRUE ~ EXENDT
  )) %>%
  ## Derive TRTDURD ----
  derive_vars_duration(TRTDURD, start_date = TRTSDT, end_date = TRTEDT) %>%
  ## Derive completion flags ----
  mutate(COMP0FL=case_when(is.na(TRTSDT) | TRTSDT>RFENDT ~ "N", T ~ "Y")) %>%
  mutate(COMP2FL=case_when(is.na(VISIT4DT) | VISIT4DT>RFENDT ~ "N", T ~ "Y")) %>%
  mutate(COMP8FL=case_when(is.na(VISIT8DT) | VISIT8DT>RFENDT ~ "N", T ~ "Y")) %>%
  mutate(COMP16FL=case_when(is.na(VISIT10DT) | VISIT10DT>RFENDT ~ "N", T ~ "Y")) %>%
  mutate(COMP24FL=case_when(is.na(VISIT12DT) | VISIT12DT>RFENDT ~ "N", T ~ "Y")) %>%
  ## Derive CUMDOSE ----
  mutate(CUMDOSE=round(case_when(
    ARMN %in% c(0, 54) & !is.na(TRTDURD) ~ as.integer(TRT01PN*TRTDURD),
    ARMN == 81 & VISIT12DT<=TRTEDT & !is.na(TRTDURD) ~ as.integer(
      54*(VISIT4DT-TRTSDT+1) +
      81*(VISIT12DT-VISIT4DT) +
      54*(TRTEDT-VISIT12DT)
    ),
    ARMN == 81 & VISIT4DT<=TRTEDT & !is.na(TRTDURD) ~ as.integer(
      54*(VISIT4DT-TRTSDT+1) +
      81*(TRTEDT-VISIT4DT)
    ),
    ARMN == 81 & TRTSDT<=TRTEDT & !is.na(TRTDURD) ~ as.integer(
      54*(TRTEDT-TRTSDT+1)
    ),
    TRUE ~ NA_integer_
  ), 1)) %>%
  ## Derive AVGDD ----
  mutate(AVGDD=round(CUMDOSE/TRTDURD, 1))

adsl2 <- adsl1 %>%
  ## Derive WEIGHTBL ----
  derive_vars_merged(
    dataset_add = vs, by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(WEIGHTBL0 = VSSTRESN),
    filter_add = VSTESTCD == "WEIGHT" & VISITNUM == 3
  ) %>%
  mutate(WEIGHTBL=round(WEIGHTBL0, 1)) %>%
  ## Derive HEIGHTBL ----
  derive_vars_merged(
    dataset_add = vs, by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(HEIGHTBL0 = VSSTRESN),
    filter_add = VSTESTCD == "HEIGHT" & VISITNUM == 1
  ) %>%
  mutate(HEIGHTBL=round(HEIGHTBL0, 1)) %>%
  ## Derive BMIBL ----
  mutate(BMIBL=round(WEIGHTBL / ((HEIGHTBL/100)**2), 1)) %>%
  ## Derive BMIBLGR1 ----
  mutate(BMIBLGR1=format_bmicat(BMIBL)) %>%
  ## Derive EDUCLVL ----
  derive_vars_merged(
    dataset_add = sc, by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(EDUCLVL = SCSTRESN),
    filter_add = SCTESTCD == "EDLEVEL"
  ) %>%
  ## Derive DISONSTDTC ----
  derive_vars_merged(
    dataset_add = mh, by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(DISONSTDTC = MHSTDTC),
    filter_add = MHCAT == "PRIMARY DIAGNOSIS"
  )

adsl3 <- adsl2 %>%
  ## Derive VISNUMEN ----
  derive_vars_merged(
    dataset_add = ds, by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(VISNMEN1 = VISITNUM),
    filter_add = DSTERM == "PROTOCOL COMPLETED"
  ) %>%
  derive_vars_merged(
    dataset_add = ds, by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(VISNMEN2 = VISITNUM),
    filter_add = DSTERM == "ADVERSE EVENT"
  ) %>%
  mutate(VISNUMEN=case_when(
    VISNMEN1 == 13 ~ 12,
    VISNMEN2 == 13 ~ 12,
    T ~ VISNMEN1
  )) %>%
  ## Derive EOSSTT ----
  mutate(EOSSTT=format_eosstt(DCDECOD)) %>%
  ## Derive DISONSDT ----
  derive_vars_dt(
    dtc = DISONSTDTC,
    new_vars_prefix = "DISONS",
    flag_imputation = "none"
  ) %>%
  ## Derive DURDIS ----
  derive_vars_duration(
    DURDIS0, start_date = DISONSDT, end_date = VISIT1DT,
    out_unit = "months"
  ) %>%
  mutate(DURDIS=round(DURDIS0, 1)) %>%
  ## Derive DURDSGR1 ----
  mutate(DURDSGR1 = format_durdisc(DURDIS))

adsl4 <- adsl3 %>%
  ## Derive EFFFL ----
  derive_var_merged_exist_flag(
    dataset_add = qs,
    by_vars = vars(STUDYID, USUBJID),
    new_var = EFFFL1,
    condition = VISITNUM > 3 & QSCAT == "ALZHEIMER'S DISEASE ASSESSMENT SCALE",
    false_value = "N",
    missing_value = "N"
  ) %>%
  derive_var_merged_exist_flag(
    dataset_add = qs,
    by_vars = vars(STUDYID, USUBJID),
    new_var = EFFFL2,
    condition = VISITNUM > 3 & QSCAT == "CLINICIAN'S INTERVIEW-BASED IMPRESSION OF CHANGE (CIBIC+)",
    false_value = "N",
    missing_value = "N"
  ) %>%
  mutate(EFFFL = case_when(EFFFL1=="Y" & EFFFL2=="Y" ~ "Y", T ~ "N")) %>%
  ## Derive SAFFL ----
  mutate(SAFFL = case_when(ITTFL=="Y" & !is.na(TRTSDT) ~ "Y", T ~ "N")) %>%
  ## Final select ----
  select(STUDYID, USUBJID, SUBJID, SITEID, SITEGR1, ARM, TRT01P, TRT01PN, TRT01A,
         TRT01AN, TRTSDT, TRTEDT, TRTDURD, AVGDD, CUMDOSE, AGE, AGEGR1, AGEGR1N,
         AGEU, RACE, RACEN, SEX, ETHNIC, SAFFL, ITTFL, EFFFL, COMP8FL, COMP16FL,
         COMP24FL, DISCONFL, DSRAEFL, DTHFL, BMIBL, BMIBLGR1, HEIGHTBL, WEIGHTBL,
         EDUCLVL, DISONSDT, DURDIS, DURDSGR1, VISIT1DT, RFSTDTC, RFENDTC, VISNUMEN,
         RFENDT, DCDECOD, EOSSTT, DCSREAS, MMSETOT)

# Add labels

adsl_labels <- c(
  'Study Identifier', 'Unique Subject Identifier', 'Subject Identifier for the Study', 'Study Site Identifier',
  'Pooled Site Group 1', 'Description of Planned Arm', 'Planned Treatment for Period 01',
  'Planned Treatment for Period 01 (N)', 'Actual Treatment for Period 01', 'Actual Treatment for Period 01 (N)',
  'Date of First Exposure to Treatment', 'Date of Last Exposure to Treatment', 'Total Treatment Duration (Days)',
  'Avg Daily Dose (as planned)', 'Cumulative Dose (as planned)', 'Age', 'Pooled Age Group 1', 'Pooled Age Group 1 (N)',
  'Age Units', 'Race', 'Race (N)', 'Sex', 'Ethnicity', 'Safety Population Flag', 'Intent-To-Treat Population Flag',
  'Efficacy Population Flag', 'Completers of Week 8 Population Flag', 'Completers of Week 16 Population Flag',
  'Completers of Week 24 Population Flag', 'Did the Subject Discontinue the Study?', 'Discontinued due to AE?',
  'Subject Died?', 'Baseline BMI (kg/m^2)', 'Pooled Baseline BMI Group 1', 'Baseline Height (cm)', 'Baseline Weight (kg)',
  'Years of Education', 'Date of Onset of Disease', 'Duration of Disease (Months)', 'Pooled Disease Duration Group 1',
  'Date of Visit 1', 'Subject Reference Start Date/Time', 'Subject Reference End Date/Time',
  'End of Trt Visit (Vis 12 or Early Term.)', 'Date of Discontinuation/Completion', 'Standardized Disposition Term',
  'End of Study Status', 'Reason for Discontinuation from Study', 'MMSE Total'
)

adsl5 <- data.frame(mapply(
  function (var, lab) {
    col <- adsl4[[var]]
    attr(col, "label") <- lab
    return(col)
  },
  var=colnames(adsl4), lab=adsl_labels,
  SIMPLIFY = FALSE
))

# Save output ----

xportr_write(adsl5, "adam/adsl.xpt")
