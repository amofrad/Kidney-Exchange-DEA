library(dplyr)
library(tidyverse)  
library(lubridate)


kidney_data <- read_delim("KIDPAN_DATA.DAT", 
                   delim = "\t", escape_double = FALSE,
                   col_names = FALSE, trim_ws = TRUE)

colnames(kidney_data) <- c("WL_ORG", "COD_WL", "COD_OSTXT_WL", "NUM_PREV_TX", "CURRENT_PRA", "PEAK_PRA", "USE_WHICH_PRA", "CREAT_CLEAR", "GFR", "DONATION", "ON_DIALYSIS", "MAX_KDPI_LOCAL_ZERO_ABDR", "MAX_KDPI_LOCAL_NON_ZERO_ABDR", "MAX_KDPI_IMPORT_ZERO_ABDR", "MAX_KDPI_IMPORT_NON_ZERO_ABDR", "C_PEPTIDE", "C_PEPTIDEDATE", "A2A2B_ELIGIBILITY", "A1", "A2", "B1", "B2", "DR1", "DR2", "ANTIBODY_TESTED", "GENDER", "ABO", "WGT_KG_TCR", "HGT_CM_TCR", "BMI_TCR", "CITIZENSHIP", "CITIZEN_COUNTRY", "PERM_STATE", "EDUCATION", "FUNC_STAT_TCR", "DGN_TCR", "DGN_OSTXT_TCR", "DGN2_TCR", "DGN2_OSTXT_TCR", "DIAB", "DRUGTRT_COPD", "TOT_SERUM_ALBUM", "C_PEPTIDE_PA_TCR", "HBA1C_PA_TCR", "SSDMF_DEATH_DATE", "INIT_CURRENT_PRA", "INIT_PEAK_PRA", "INIT_STAT", "INIT_WGT_KG", "INIT_HGT_CM", "INIT_CPRA", "END_CPRA", "INIT_EPTS", "END_EPTS", "REM_CD", "DAYSWAIT_CHRON", "END_STAT", "INIT_AGE","ACTIVATE_DATE", "CREAT_CLEAR_DATE", "DEATH_DATE", "DIALYSIS_DATE", "END_DATE", "GFR_DATE", "INIT_DATE", "WT_QUAL_DATE", "ETHNICITY", "ETHCAT", "PT_CODE","INIT_BMI_CALC", "END_BMI_CALC", "DAYSWAIT_ALLOC", "COMPOSITE_DEATH_DATE", "WLHR", "WLHL", "WLIN", "WLKI", "WLKP", "WLLI", "WLLU", "WLPA", "WLPI", "WLVC", "REGION", "INACT_REASON_CD", "BW4", "BW6", "C1", "C2", "DR51", "DR51_2", "DR52", "DR52_2", "DR53", "DR53_2", "DQ1", "DQ2", "WL_ID_CODE", "PERIP_VASC", "EXH_PERIT_ACCESS", "AGE_DIAB", "EXH_VASC_ACCESS", "YR_ENTRY_US_TCR", "WORK_INCOME_TCR", "ACADEMIC_PRG_TCR", "ACADEMIC_LEVEL_TCR", "MALIG_TCR_KI", "PRI_PAYMENT_TCR_KI", "MALIG_TCR_PA", "PRI_PAYMENT_TCR_PA", "PREV_TX", "PREV_KI_TX", "PREV_PA_TX", "ACADEMIC_LEVEL_TRR", "ACADEMIC_PRG_TRR", "FUNC_STAT_TRR", "MALIG_TRR", "MALIG_OSTXT_TRR", "MALIG_TY_TRR", "PERM_STATE_TRR", "PRI_PAYMENT_TRR_KI", "PRI_PAYMENT_CTRY_TRR_KI", "WORK_INCOME_TRR", "TX_DATE", "ACUTE_REJ_EPI_KI", "CREAT_TRR", "FIRST_WK_DIAL", "ORG_REC_ON", "PREV_PREG", "REC_ON_ICE", "REC_ON_PUMP", "SERUM_CREAT", "L_FIN_FLOW_RATE_TX", "L_FIN_RESIST_TX", "R_FIN_FLOW_RATE_TX", "R_FIN_RESIST_TX", "PRE_TX_TXFUS", "FIN_RESIST_TX", "TXHRT", "TXINT", "TXKID", "TXLIV", "TXLNG", "TXPAN", "TXVCA", "PRI_PAYMENT_CTRY_TCR_KI", "PREV_MALIG_TY", "PREV_MALIG_TY_OSTXT", "RDA1", "RDA2", "RDB1", "RDB2", "RDDR1", "RDDR2", "DON_RETYP", "RESUM_MAINT_DIAL_DT", "DA1", "DA2", "DB1", "DB2", "DDR1", "DDR2", "RA1", "RA2", "RB1", "RB2", "RDR1", "RDR2", "AMIS", "BMIS", "DRMIS", "HLAMIS", "NPKID", "NPPAN", "END_CPRA_DETAIL", "HAPLO_TY_MATCH_DON", "AGE_DON", "DDAVP_DON", "CMV_OLD_LIV_DON", "CMV_DON", "CMV_TEST_DON", "EBV_TEST_DON", "HBV_TEST_DON", "HCV_TEST_DON", "CMV_NUCLEIC_DON", "CMV_IGG_DON", "CMV_IGM_DON", "EBV_DNA_DON", "EBV_IGG_DON", "EBV_IGM_DON", "HBV_CORE_DON", "HBV_SUR_ANTIGEN_DON", "ETHCAT_DON", "COD_CAD_DON", "DEATH_CIRCUM_DON", "DEATH_MECH_DON", "CITIZENSHIP_DON", "HEP_C_ANTI_DON", "HCV_RNA_DON", "ABO_DON", "DON_TY", "GENDER_DON", "HOME_STATE_DON", "WARM_ISCH_TM_DON", "HCV_RIBA_DON", "HCV_ANTIBODY_DON", "LIV_DON_TY", "CITIZEN_COUNTRY_DON", "COD_OSTXT_DON", "CONTROLLED_DON", "CORE_COOL_DON", "NON_HRT_DON", "ANTIHYPE_DON", "BLOOD_INF_DON", "BLOOD_INF_CONF_DON", "BUN_DON", "CREAT_DON", "DOBUT_DON_OLD", "DOPAMINE_DON_OLD", "HTLV1_OLD_DON", "HTLV2_OLD_DON", "OTH_DON_MED1_OSTXT_DON_OLD", "OTH_DON_MED2_OSTXT_DON_OLD", "OTH_DON_MED3_OSTXT_DON_OLD", "OTHER_INF_DON", "OTHER_INF_CONF_DON", "OTHER_INF_OSTXT_DON", "PRETREAT_MED_DON_OLD", "PT_DIURETICS_DON", "PT_STEROIDS_DON", "PT_T3_DON", "PT_T4_DON", "PT_OTH2_OSTXT_DON", "PT_OTH3_OSTXT_DON", "PT_OTH4_OSTXT_DON", "PT_OTH1_OSTXT_DON", "PULM_INF_DON", "PULM_INF_CONF_DON", "SGOT_DON", "SGPT_DON", "TBILI_DON", "URINE_INF_DON", "URINE_INF_CONF_DON", "VASODIL_DON", "VDRL_DON", "CLIN_INFECT_DON", "HYPERTENS_DUR_DON", "CANCER_FREE_INT_DON", "CANCER_OTH_OSTXT_DON", "CONTIN_ALCOHOL_OLD_DON", "CONTIN_CIG_DON", "CONTIN_IV_DRUG_OLD_DON", "CONTIN_COCAINE_DON", "CONTIN_OTH_DRUG_DON", "DIET_DON", "DIURETICS_DON", "EXTRACRANIAL_CANCER_DON"         , "HIST_ALCOHOL_OLD_DON", "CANCER_SITE_DON", "HIST_CIG_DON", "DIABDUR_DON", "HIST_COCAINE_DON", "HIST_HYPERTENS_DON", "HIST_IV_DRUG_OLD_DON", "INSULIN_DEP_DON", 
                    "INTRACRANIAL_CANCER_DON", "OTHER_HYPERTENS_MED_DON", "HIST_CANCER_DON", "HIST_INSULIN_DEP_DON", "INSULIN_DUR_DON", "HIST_DIABETES_DON", "HIST_OTH_DRUG_DON", "SKIN_CANCER_DON", "DIABETES_DON", "LIV_DON_TY_OSTXT", "HEPARIN_DON", "ARGININE_DON", "INSULIN_DON", "HGT_CM_DON_CALC", "WGT_KG_DON_CALC", "BMI_DON_CALC", "KDPI", "KDRI_MED", "KDRI_RAO", "HBV_NAT_DON", "HCV_NAT_DON", "HIV_NAT_DON", "END_STAT_KI", "CREAT6M", "CREAT1Y", "DIAL_DATE", "RETXDATE_KI", "FAILDATE_KI", "PUMP_KI", "ABO_MAT", "AGE", "DISTANCE", "RESUM_MAINT_DIAL", "DIAL_TRR", "DIAG_KI", "DIAG_OSTXT_KI", "COLD_ISCH_KI", "GRF_STAT_KI", "GRF_FAIL_CAUSE_OSTXT_KI", "GRF_FAIL_CAUSE_TY_KI", "DWFG_KI", "PRVTXDIF_KI", "GTIME_KI", "GSTATUS_KI", "COD_KI", "COD_OSTXT_KI", "COD2_KI", "COD2_OSTXT_KI", "COD3_KI", "COD3_OSTXT_KI", "DAYSWAIT_CHRON_KI", "TX_PROCEDUR_TY_KI", "TRTREJ1Y_KI", "TRTREJ6M_KI", "MULTIORG", "PRI_PAYMENT_TRR_PA", "PRI_PAYMENT_CTRY_TRR_PA", "ART_RECON", "ART_RECON_OSTXT", "DUCT_MGMT", "DUCT_MGMT_OSTXT", "GRF_PLACEM", "PRE_AVG_INSULIN_USED_TRR", "PRE_AVG_INSULIN_USED_OLD_TRR", "ACUTE_REJ_EPI_PA", "PA_PRESERV_TM", "VASC_MGMT", "VEN_EXT_GRF", "INSULIN_PA", "INSULIN_RESUMED_DATE_PA", "INSULIN_DOSAGE_PA", "INSULIN_DURATION_PA", "METHOD_BLOOD_SUGAR_CONTROL_PA", "BLOOD_SUGAR_MEDICATION_PA", "BLOOD_SUGAR_MED_RESUMED_DATE_PA", "BLOOD_SUGAR_DIET_PA", "C_PEPTIDE_PA_TRR", "HBA1C_PA_TRR", "INSULIN_DOSAGE_OLD_PA", "PRI_PAYMENT_CTRY_TCR_PA", "PK_DA1", "PK_DA2", "PK_DB1", "PK_DB2", "PK_DDR1", "PK_DDR2", "ENTERIC_DRAIN", "ENTERIC_DRAIN_DT", "END_STAT_PA", "FAILDATE_PA", "DIAG_PA", "DIAG_OSTXT_PA", "GRF_STAT_PA", "GRF_FAIL_CAUSE_OSTXT_PA", "GRF_FAIL_CAUSE_TY_PA", "OTH_GRF_FAIL_CAUSE_OSTXT_PA", "GRF_VASC_THROMB_PA", "INFECT_PA", "BLEED_PA", "ANAST_LK_PA", "REJ_ACUTE_PA", "REJ_HYPER_PA", "BIOP_ISLET_PA", "PANCREATIT_PA", "REJ_CHRONIC_PA", "PX_NON_COMPL_PA", "RETXDATE_PA", "PRVTXDIF_PA", "GTIME_PA", "GSTATUS_PA", "COD_PA", "COD_OSTXT_PA", "COD2_PA", "COD2_OSTXT_PA", "COD3_PA", "COD3_OSTXT_PA", "DAYSWAIT_CHRON_PA", "TX_PROCEDUR_TY_PA", "TRTREJ1Y_PA", "TRTREJ6M_PA", "ORGAN", "CMV_IGG", "CMV_IGM", "EBV_SEROSTATUS", "HBV_CORE", "HBV_SUR_ANTIGEN", "HCV_SEROSTATUS", "HIV_SEROSTATUS", "CMV_STATUS", "HBV_SURF_TOTAL", "HIV_NAT", "HCV_NAT", "HBV_NAT", "PREV_TX_ANY", "PREV_TX_ANY_N", "TX_TYPE", "MED_COND_TRR", "PX_STAT", "PX_STAT_DATE", "PREV_KI_DATE", "FUNC_STAT_TRF", "SHARE_TY", "PSTATUS", "PTIME", "LOS", "PAYBACK", "ECD_DONOR", "AGE_GROUP", "MALIG", "MALIG_TY_OSTXT", "MALIG_TY", "HGT_CM_CALC", "WGT_KG_CALC", "BMI_CALC", "STATUS_TCR", "STATUS_TRR", "STATUS_DDR", "VAL_DT_DDR", "STATUS_LDR", "VAL_DT_LDR", "VAL_DT_TCR", "VAL_DT_TRR", "LT_ONE_WEEK_DON", "REJ_BIOPSY", "REJCNF_KI", "REJTRT_KI", "REJCNF_PA", "REJTRT_PA", "TRR_ID_CODE", "ADMISSION_DATE", "DISCHARGE_DATE", "COMPL_ABSC", "COMPL_ANASLK", "COMPL_PANCREA", "OTH_COMPL_OSTXT", "SURG_INCIS", "OPER_TECH", "EDUCATION_DON", "KI_CREAT_PREOP", "KI_PROC_TY", "PRI_PAYMENT_DON", "PRI_PAYMENT_CTRY_DON", "MEDICARE_DON", "MEDICAID_DON", "OTH_GOVT_DON", "PRIV_INS_DON", "HMO_PPO_DON", "SELF_DON", "DONATION_DON", "FREE_DON", "RECOV_OUT_US", "RECOV_COUNTRY", "PROTEIN_URINE", "LIPASE", "AMYLASE", "INOTROP_AGENTS", "CARDARREST_NEURO", "RESUSCIT_DUR", "INOTROP_SUPPORT_DON", "TATTOOS", "LT_KI_BIOPSY", "LT_KI_GLOMERUL", "RT_KI_BIOPSY", "RT_KI_GLOMERUL", "REFERRAL_DATE", "RECOVERY_DATE", "ADMIT_DATE_DON", "DONOR_ID", "HBSAB_DON", "EBV_IGG_CAD_DON", "EBV_IGM_CAD_DON", "HBV_DNA_DON", "CDC_RISK_HIV_DON", "INO_PROCURE_AGENT_1", "INO_PROCURE_AGENT_2", "INO_PROCURE_AGENT_3", "INO_PROCURE_OSTXT_1", "INO_PROCURE_OSTXT_2", "INO_PROCURE_OSTXT_3", "DATA_TRANSPLANT", "DATA_WAITLIST", "CTR_CODE", "OPO_CTR_CODE", "INIT_OPO_CTR_CODE", "END_OPO_CTR_CODE", "LISTING_CTR_CODE")


################################################################################
################################# DATA PREP ####################################
################################################################################

kidney_data <- kidney_data %>%
  # KPD only
  filter(LIV_DON_TY == 9) %>%
  mutate(
    # Transplant date and year
    TX_YEAR = year(mdy(TX_DATE)),
    
    # Clean waitlist duration (".") → NA → numeric
    DAYSWAIT_CHRON_KI = ifelse(DAYSWAIT_CHRON_KI == ".", 
                               NA, 
                               DAYSWAIT_CHRON_KI),
    DAYSWAIT_CHRON_KI = as.numeric(DAYSWAIT_CHRON_KI),
    
    # Dialysis time (same logic as before, just inside mutate)
    DIAL_TIME = ifelse(
      DIALYSIS_DATE == ".",
      0,
      as.numeric(
        as.Date(INIT_DATE, "%m/%d/%Y") -
          as.Date(DIALYSIS_DATE, "%m/%d/%Y")
      )
    ),
    
    # Recipient ethnicity recode
    ETHCAT = case_when(
      ETHCAT == 1   ~ "White",
      ETHCAT == 2   ~ "Black",
      ETHCAT == 4   ~ "Hispanic",
      ETHCAT == 5   ~ "Asian",
      ETHCAT %in% c(6, 7, 9, 998) ~ "Other",
      TRUE          ~ as.character(ETHCAT)
    ),
    
    # Donor ethnicity recode
    ETHCAT_DON = case_when(
      ETHCAT_DON == 1   ~ "White",
      ETHCAT_DON == 2   ~ "Black",
      ETHCAT_DON == 4   ~ "Hispanic",
      ETHCAT_DON == 5   ~ "Asian",
      ETHCAT_DON %in% c(6, 7, 9, 998) ~ "Other",
      TRUE              ~ as.character(ETHCAT_DON)
    ),
    
    # Ethnicity indicators
    ETH_White    = as.integer(ETHCAT == "White"),
    ETH_Asian    = as.integer(ETHCAT == "Asian"),
    ETH_Black    = as.integer(ETHCAT == "Black"),
    ETH_Hispanic = as.integer(ETHCAT == "Hispanic")
  ) %>%
  # Restrict to analysis years
  filter(TX_YEAR >= 2010, TX_YEAR <= 2019) %>%
  # Factor-ise ethnicity variables
  mutate(
    ETHCAT     = as.factor(ETHCAT),
    ETHCAT_DON = as.factor(ETHCAT_DON)
  )


Donors <- read_delim("LIVING_DONOR_DATA.DAT", 
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE)

colnames(Donors) <- c("REOP_BILIARY", "REOP_BILIARY_DT", "DON_DATE", "AGE_DON", "ETHCAT_DON", "REGION", "LUNG_RECOV", "KIDNEY_RECOV", "LIVER_RECOV", "DON_ORG2", "CITIZENSHIP", "LIV_DON_TY", "LIV_DON_TY_OSTXT", "EDUCATION", "KI_CREAT_PREOP", "BP_PREOP_SYST", "BP_PREOP_DIAST", "COD", "KI_PROC_TY", "MARITAL_STAT", "HEALTH_INS", "FUNC_STAT", "PHYSICAL_CAPACITY", "WORK_INCOME", "HIST_CANCER", "HIST_CANCER_OSTXT", "CANCER_FREE", "COD_OSTXT", "HIST_HYPER", "HYPER_DIET", "HYPER_DIUR", "HYPER_MEDS", "PREOP_URINE_RATIO", "DIABETES", "MACRO_FAT", "MICRO_FAT", "CONVERT_OPEN_KI", "NON_AUTO_BLOOD", "PRBC_UNITS", "PLATELETS_UNITS", "FFP_UNITS", "VASC_COMP_KI", "VASC_COMP_KI_INTER", "VASC_COMP_KI_INTER_OSTXT", "OTH_COMP_KI", "OTH_COMP_KI_INTER", "OTH_COMP_KI_INTER_OSTXT", "REOPERATION_KI", "REOP_BLEED_KI", "REOP_HERNIA_KI", "REOP_BOWEL_KI", "REOP_VASC_KI", "REOP_OTH_KI", "REOP_OTH_KI_OSTXT", "READMISSION_KI", "READMISSION_KI_REASON", "READMISSION_KI_OSTXT", "OTH_INTER_PROC_KI", "OTH_INTER_PROC_KI_OSTXT", "KI_CREAT_POSTOP", "BP_POSTOP_SYST", "BP_POSTOP_DIAST", "HYPERTENSION", "POSTOP_URINE_RATIO", "PREOP_BILI", "PREOP_SGOT_AST", "PREOP_ALK_PHOS", "PREOP_ALBUM", "PREOP_INR", "BIOPSY_LI", "LI_PROC_TY", "BILIARY_COMP", "BILIARY_COMP_GRADE", "VASC_COMP_LI", "VASC_COMP_LI_INTER", "VASC_COMP_LI_INTER_OSTXT", "OTH_COMP_LI", "OTH_COMP_LI_INTER", "OTH_COMP_LI_INTER_OSTXT", "REOPERATION_LI", "REOP_BLEED_LI", "REOP_HERNIA_LI", "REOP_BOWEL_LI", "REOP_VASC_LI", "REOP_OTH_LI", "REOP_OTH_LI_OSTXT", "REOP_LI_FAIL", "READMISSION_LI", "READMISSION_LI_REASON", "READMISSION_LI_OSTXT", "OTH_INTER_PROC_LI", "OTH_INTER_PROC_LI_OSTXT", "POSTOP_SGOT_AST", "POSTOP_ALK_PHOS", "POSTOP_ALBUM", "POSTOP_CREAT_LI", "POSTOP_INR", "POSTOP_SGPT_ALT", "POSTOP_BILI", "PREOP_FVC_BEFORE", "PREOP_FVC_AFTER", "PREOP_FEV1_BEFORE", "PREOP_FEV1_AFTER", "PREOP_FEF_BEFORE", "PREOP_FEF_AFTER", "PREOP_TLC_BEFORE", "PREOP_TLC_AFTER", "PREOP_LUNG_CAP", "PREOP_PAO2", "HIST_CIG", "PACK_YRS", "DUR_ABSTINENCE", "TOBACCO_USE", "LU_PROC_TY", "CONVERT_OPEN_LU", "INTRAOP_COMP", "INTRAOP_COMP_REASON", "SACRIFICE_LOBE", "ARRHYTHMIA", "ANESTHETIC_COMP", "INTRAOP_COMP_OSTXT", "LU_COMP", "LU_COMP_REASON", "LU_COMP_OSTXT", "THORAC_TUBES", "ARRHYTHMIA_POSTOP", "READMISSION_LU", "READMISSION_LU_REASON", "READMISSION_LU_OSTXT", "PREDON_HGT", "PREDON_WGT", "PREOP_SGPT_ALT", "PREOP_CREAT_LI", "PREOP_URINE_PROTEIN", "POSTOP_URINE_PROTEIN", "PX_STAT", "CITIZEN_COUNTRY", "DEATH_DT", "ORG_RECOVERY_DT", "INIT_DISCHARGE_DT", "REOP_BLEED_KI_DT", "REOP_HERNIA_KI_DT", "REOP_BOWEL_KI_DT", "REOP_VASC_KI_DT", "REOP_OTH_KI_DT", "READMISSION_KI_DT", "OTH_INTER_PROC_KI_DT", "POSTOP_TEST_DT", "REOP_BLEED_LI_DT", "REOP_HERNIA_LI_DT", "REOP_BOWEL_LI_DT", "REOP_VASC_LI_DT", "REOP_OTH_LI_DT", "REOP_LI_FAIL_DT", "READMISSION_LI_DT", "OTH_INTER_PROC_LI_DT", "READMISSION_LU_DT", "CMV_IGG", "CMV_IGM", "CMV_NUCLEIC", "EBV_IGG", "EBV_IGM", "HBV_CORE", "HBV_DNA", "HBV_SUR_ANTIGEN", "HCV_ANTIBODY", "HCV_RIBA", "HCV_RNA", "VIRUSES_TESTED", "CMV_TOTAL", "EBV_TOTAL", "HOME_STATE", "GENDER", "ABO", "YR_ENTRY_US", "WGT_KG", "DON_ORG", "STATUS_LDR", "VAL_DT_LDR", "DBW4", "DBW6", "DC1", "DC2", "DDP1", "DDP2", "DDPA1", "DDPA2", "DDR51", "DDR51_2", "DDR52", "DDR52_2", "DDR53", "DDR53_2", "DDQ1", "DDQ2", "DDQA1", "DDQA2", "RECOV_FACILITY_CODE", "DONOR_ID")
## Store the original kidney_data variable names before the join
base_vars <- names(kidney_data)

## Vectorized helpers (as before)
calculate_eGFR_vec <- function(KI_CREAT_PREOP, AGE_DON, GENDER_DON, African_American) {
  kappa <- ifelse(GENDER_DON == "F", 0.7, 0.9)
  alpha <- ifelse(GENDER_DON == "F", -0.329, -0.411)
  
  base <- 141 *
    pmin(KI_CREAT_PREOP / kappa, 1)^alpha *
    pmax(KI_CREAT_PREOP / kappa, 1)^(-1.209) *
    0.993^AGE_DON
  
  base <- ifelse(GENDER_DON == "F", base * 1.018, base)
  base <- ifelse(African_American, base * 1.159, base)
  base
}

calculate_LKDPI_vec <- function(age, eGFR, BMI, african_american,
                                history_of_cigarette_use, SBP, both_male,
                                ABO_incompatible, unrelated,
                                HLA_B_mismatches, HLA_DR_mismatches, D_RWR) {
  LKDPI <- rep(-11.30, length(age))
  
  LKDPI <- LKDPI + ifelse(!is.na(age) & age > 50, 1.85 * (age - 50), 0)
  LKDPI <- LKDPI - 0.381 * eGFR + 1.17 * BMI
  LKDPI <- LKDPI + ifelse(african_american,        22.34, 0)
  LKDPI <- LKDPI + ifelse(history_of_cigarette_use, 14.33, 0)
  LKDPI <- LKDPI + 0.44 * SBP
  LKDPI <- LKDPI + ifelse(both_male,              -21.68, 0)
  LKDPI <- LKDPI + ifelse(ABO_incompatible,        27.30, 0)
  LKDPI <- LKDPI + ifelse(unrelated,              -10.61, 0)
  LKDPI <- LKDPI + 8.57 * HLA_B_mismatches + 8.26 * HLA_DR_mismatches
  LKDPI <- LKDPI - 50.87 * pmin(D_RWR, 0.9)
  
  LKDPI
}

## Clean Donors, join, compute eGFR + LKDPI, trim columns, and recode types
dmu_data <- kidney_data %>%
  mutate(DONOR_ID = as.character(DONOR_ID)) %>%
  # join donor info
  left_join(
    Donors %>%
      mutate(
        DONOR_ID        = as.character(DONOR_ID),  # <- key fix here
        KI_CREAT_PREOP   = as.numeric(KI_CREAT_PREOP),
        African_American = ETHCAT_DON == 2
      ),
    by = "DONOR_ID"
  ) %>%
  arrange(DONOR_ID) %>%
  # compute eGFR + LKDPI ingredients + LKDPI itself
  mutate(
    age = as.numeric(AGE_DON.x),
    
    eGFR = ifelse(
      complete.cases(KI_CREAT_PREOP.y, AGE_DON.y, GENDER.y, African_American),
      calculate_eGFR_vec(
        KI_CREAT_PREOP.y,
        AGE_DON.y,
        GENDER.y,
        African_American
      ),
      NA_real_
    ),
    
    BMI                      = as.numeric(BMI_DON_CALC),
    african_american         = (ETHCAT_DON.x == "Black"),
    history_of_cigarette_use = (HIST_CIG_DON == "Y"),
    SBP                      = as.numeric(BP_PREOP_SYST),
    both_male                = (GENDER.x == "M" & GENDER_DON == "M"),
    ABO_incompatible         = (ABO_MAT == 3),
    unrelated                = !(LIV_DON_TY.x %in% 1:6),
    HLA_B_mismatches         = as.integer(B1 == DB1)  + as.integer(B2 == DB2),
    HLA_DR_mismatches        = as.integer(DDR1 == DR1) + as.integer(DDR2 == DR2),
    D_RWR                    = as.numeric(WGT_KG_DON_CALC) / as.numeric(WGT_KG_CALC),
    
    LKDPI = calculate_LKDPI_vec(
      age,
      eGFR,
      BMI,
      african_american,
      history_of_cigarette_use,
      SBP,
      both_male,
      ABO_incompatible,
      unrelated,
      HLA_B_mismatches,
      HLA_DR_mismatches,
      D_RWR
    )
  ) %>%
  # keep only rows with non-missing LKDPI
  filter(!is.na(LKDPI)) %>%
  # drop ".x" suffix from original kidney variables
  rename_with(~ sub("\\.x$", "", .x)) %>%
  # keep original kidney variables + new LKDPI-related ones
  dplyr::select(any_of(c(
    base_vars,
    "LKDPI", "eGFR", "SBP", "HLA_B_mismatches", "HLA_DR_mismatches", "ABO_incompatible"
  ))) %>%
  # ON_DIALYSIS and core numeric/factor casting
  filter(!is.na(ON_DIALYSIS)) %>%
  mutate(
    ON_DIALYSIS   = factor(ifelse(ON_DIALYSIS == "Y", 1L, 0L)),
    GTIME_KI      = as.numeric(GTIME_KI),
    ABO           = factor(ABO),
    AGE           = as.numeric(AGE),
    GENDER        = factor(GENDER),
    PRA           = as.numeric(END_CPRA),
    WGT_KG_CALC   = as.numeric(WGT_KG_CALC),
    BMI           = as.numeric(BMI_CALC),
    DIABETES_DON  = factor(DIABETES_DON)
  ) %>%
  # restrict to the 4 focal ethnic groups
  filter(ETHCAT %in% c("Asian", "Black", "Hispanic", "White")) %>%
  mutate(
    ETHCAT     = droplevels(ETHCAT),
    ETHCAT_DON = droplevels(ETHCAT_DON)
  ) %>%
  # recode PRI_PAYMENT_TRR_KI, ACUTE_REJ_EPI_KI, EDUCATION
  mutate(
    PRI_PAYMENT_TRR_KI = case_when(
      PRI_PAYMENT_TRR_KI %in% 1:14 ~ as.factor(PRI_PAYMENT_TRR_KI),
      TRUE ~ as.factor("Other")
    ),
    ACUTE_REJ_EPI_KI = case_when(
      ACUTE_REJ_EPI_KI %in% 1:3 ~ as.factor(ACUTE_REJ_EPI_KI),
      TRUE ~ as.factor("Unknown")
    ),
    EDUCATION = case_when(
      EDUCATION %in% c(1:6, 996, 998) ~ as.factor(EDUCATION),
      TRUE ~ as.factor("Unknown")
    )
  ) %>%
  mutate(
    PRI_PAYMENT_TRR_KI = droplevels(PRI_PAYMENT_TRR_KI),
    ACUTE_REJ_EPI_KI   = droplevels(ACUTE_REJ_EPI_KI),
    EDUCATION          = droplevels(EDUCATION),
    
    AGE_DON           = as.numeric(AGE_DON),
    HGT_CM_DON_CALC   = as.numeric(HGT_CM_DON_CALC),
    WGT_KG_DON_CALC   = as.numeric(WGT_KG_DON_CALC),
    DISTANCE          = as.numeric(DISTANCE),
    LOS               = as.numeric(LOS),
    HGT_CM_CALC       = as.numeric(HGT_CM_CALC),
    
    HIST_HYPERTENS_DON = as.factor(HIST_HYPERTENS_DON),
    PREV_KI_TX         = as.factor(PREV_KI_TX),
    HIST_CANCER_DON    = as.factor(HIST_CANCER_DON),
    HIST_CIG_DON       = as.factor(HIST_CIG_DON),
    DIAG_KI            = as.factor(DIAG_KI),
    ACADEMIC_LEVEL_TCR = as.factor(ACADEMIC_LEVEL_TCR),
    DIAB               = as.factor(DIAB),
    EDUCATION          = as.factor(EDUCATION),
    ABO_MAT            = as.factor(ABO_MAT),
    PREV_TX_ANY        = as.factor(PREV_TX_ANY),
    REGION             = as.factor(REGION),
    WORK_INCOME_TCR    = as.factor(WORK_INCOME_TCR),
    MED_COND_TRR       = as.factor(MED_COND_TRR),
    TX_PROCEDUR_TY_KI  = as.factor(TX_PROCEDUR_TY_KI),
    INIT_STAT          = as.factor(INIT_STAT),
    END_STAT           = as.factor(END_STAT),
    
    B1              = as.factor(B1),
    DB1             = as.factor(DB1),
    B2              = as.factor(B2),
    DB2             = as.factor(DB2),
    DDR1            = as.factor(DDR1),
    DR1             = as.factor(DR1),
    DR2             = as.factor(DR2),
    DDR2            = as.factor(DDR2),
    ABO_incompatible = as.factor(ABO_incompatible)
  ) %>%
  # ABO dummies + GRF_FAIL_CAUSE_TY_KI + DWFG_KI + payments, citizenship, COD_KI, NA-handling, etc.
  filter(ABO %in% c("A", "B", "AB", "O")) %>%
  mutate(
    ABO = droplevels(as.factor(ABO)),
    
    # ABO indicator variables
    ABO_O  = as.numeric(ABO == "O"),
    ABO_AB = as.numeric(ABO == "AB"),
    ABO_B  = as.numeric(ABO == "B"),
    ABO_A  = as.numeric(ABO == "A"),
    
    
    # GRF failure cause labels
    GRF_FAIL_CAUSE_TY_KI = case_when(
      GRF_FAIL_CAUSE_TY_KI == 1   ~ "Hyperacute Rejection",
      GRF_FAIL_CAUSE_TY_KI == 2   ~ "Acute Rejection",
      GRF_FAIL_CAUSE_TY_KI == 3   ~ "Primary Failure",
      GRF_FAIL_CAUSE_TY_KI == 4   ~ "Graft Thrombosis",
      GRF_FAIL_CAUSE_TY_KI == 5   ~ "Infection",
      GRF_FAIL_CAUSE_TY_KI == 6   ~ "Surgical Complications",
      GRF_FAIL_CAUSE_TY_KI == 7   ~ "Urological Complications",
      GRF_FAIL_CAUSE_TY_KI == 8   ~ "Recurrent Disease",
      GRF_FAIL_CAUSE_TY_KI == 9   ~ "Primary Non-Function (Graft Never Functioned Post-Transplant)",
      GRF_FAIL_CAUSE_TY_KI == 10  ~ "Chronic Rejection",
      GRF_FAIL_CAUSE_TY_KI == 11  ~ "BK (Polyoma) Virus",
      GRF_FAIL_CAUSE_TY_KI == 12  ~ "Primary Non-Function (Graft Never Functioned Post-Transplant)",
      GRF_FAIL_CAUSE_TY_KI == 999 ~ "Other",
      TRUE                        ~ "Unknown"
    ),
    GRF_FAIL_CAUSE_TY_KI = factor(GRF_FAIL_CAUSE_TY_KI),
    
    # Death with functioning graft
    DWFG_KI = case_when(
      DWFG_KI == "Y" ~ 1L,
      DWFG_KI == "N" ~ 0L,
      TRUE           ~ NA_integer_
    ),
    DWFG_KI    = factor(DWFG_KI),
    GSTATUS_KI = as.numeric(GSTATUS_KI),
    
    BMI_DON_CALC = as.numeric(BMI_DON_CALC),
    PERM_STATE   = as.factor(PERM_STATE),
    CITIZENSHIP  = as.factor(CITIZENSHIP),
    PRI_PAYMENT_TCR_KI = as.factor(PRI_PAYMENT_TCR_KI)
  ) %>%
  mutate(
    # Primary payment type at listing
    PRI_PAYMENT_TCR_KI = case_when(
      PRI_PAYMENT_TCR_KI %in% c(3, 4, 13) ~ "Medicare",
      PRI_PAYMENT_TCR_KI == 1 ~ "Private Insurance",
      PRI_PAYMENT_TCR_KI == 2 ~ "Medicaid",
      PRI_PAYMENT_TCR_KI == 5 ~ "CHIP",
      PRI_PAYMENT_TCR_KI == 6 ~ "VA",
      PRI_PAYMENT_TCR_KI == 7 ~ "Other Gov",
      PRI_PAYMENT_TCR_KI == 8 ~ "Self",
      PRI_PAYMENT_TCR_KI %in% c(9, 10, 11, 12) ~ "Other",
      TRUE ~ NA_character_
    ),
    PRI_PAYMENT_TCR_KI = factor(PRI_PAYMENT_TCR_KI)
  ) %>%
  mutate(
    # Citizenship recode
    CITIZENSHIP = case_when(
      CITIZENSHIP == 1 ~ "US Citizen",
      CITIZENSHIP == 2 ~ "RESIDENT ALIEN",
      CITIZENSHIP %in% c(3, 5, 6) ~ "NON-RESIDENT ALIEN",
      CITIZENSHIP == 4 ~ "Non-Citizen, US Resident",
      TRUE ~ NA_character_
    ),
    CITIZENSHIP = factor(CITIZENSHIP),
    
    LISTING_CTR_CODE = as.factor(LISTING_CTR_CODE),
    COD_KI           = as.factor(COD_KI)
  ) %>%
  mutate(
    # Cause of death recode
    COD_KI = case_when(
      COD_KI %in% c("3200", "3201", "3202", "3203", "3204", "3299") ~ "Graft Fail", 
      COD_KI %in% c("3300", "3301", "3302", "3303", "3304", "3305", "3306", "3307", "3308", "3399") ~ "Infection",
      COD_KI %in% c("3400", "3401", "3402", "3499") ~ "Cardiovascular",
      COD_KI %in% c("3500", "3599") ~ "Cerebrovascular",
      COD_KI %in% c("3600", "3601", "3699") ~ "Hemorrhage",
      COD_KI %in% c("3700", "3701", "3702", "3799") ~ "Malignancy",
      COD_KI %in% c("3800", "3899") ~ "Trauma",
      COD_KI %in% c(
        "3900", "3901", "3902", "3903", "3904", "3905", "3906", "3907", "3908",
        "3909", "3910", "3911", "3912", "3913", "3914"
      ) ~ "Misc.",
      COD_KI == "3915" ~ "Primary non-function",
      COD_KI %in% c("3916", "3917") ~ "Viral Infection",
      COD_KI == "." ~ "Unknown"
    ),
    #COD_KI   = factor(COD_KI),
    CTR_CODE = as.factor(CTR_CODE)
  ) %>%
  mutate(
    TRTREJ1Y_KI   = fct_na_value_to_level(as.factor(TRTREJ1Y_KI), level = "U"),
    TRTREJ6M_KI   = fct_na_value_to_level(as.factor(TRTREJ6M_KI), level = "U"),
    WORK_INCOME_TCR = fct_na_value_to_level(as.factor(WORK_INCOME_TCR), level = "U"),
    WORK_INCOME_TRR = fct_na_value_to_level(as.factor(WORK_INCOME_TRR), level = "U"),
    FUNC_STAT_TRF   = as.factor(FUNC_STAT_TRF)
  ) %>% 
  transmute(
    Group           = ETHCAT,
    ETHCAT          = ETHCAT,
    WaitlistDuration = DAYSWAIT_CHRON_KI,
    QualityScore    = LKDPI,
    OutcomeScore    = GTIME_KI,
    
    WL_days         = DAYSWAIT_CHRON_KI,
    LKDPI           = LKDPI,
    GTIME_KI        = GTIME_KI,
    Year            = TX_YEAR,
    DIABETES_DON    = DIABETES_DON,
    eGFR            = eGFR,
    BMI_CALC        = BMI,
    WGT_KG_CALC     = WGT_KG_CALC,
    ABO             = ABO,
    ON_DIALYSIS     = ON_DIALYSIS,
    AGE             = AGE,
    AGE_DON         = AGE_DON,
    HGT_CM_DON_CALC = HGT_CM_DON_CALC,
    WGT_KG_DON_CALC = WGT_KG_DON_CALC,
    ETHCAT_DON      = ETHCAT_DON,
    HIST_HYPERTENS_DON = HIST_HYPERTENS_DON,
    PRA             = PRA,
    GENDER          = GENDER,
    DISTANCE        = DISTANCE,
    PRI_PAYMENT_TRR_KI = PRI_PAYMENT_TRR_KI,
    PREV_KI_TX      = PREV_KI_TX,
    HIST_DIABETES_DON = HIST_DIABETES_DON,
    CTR_CODE        = CTR_CODE,
    DIAG_KI         = DIAG_KI,
    BMI_DON_CALC    = BMI_DON_CALC,
    GRF_FAIL_CAUSE_TY_KI = GRF_FAIL_CAUSE_TY_KI,
    ACADEMIC_LEVEL_TCR = ACADEMIC_LEVEL_TCR,
    DIAB            = DIAB,
    EDUCATION       = EDUCATION,
    ABO_MAT         = ABO_MAT,
    PREV_TX_ANY     = PREV_TX_ANY,
    REGION          = REGION,
    MED_COND_TRR    = MED_COND_TRR,
    LOS             = LOS,
    INIT_STAT       = INIT_STAT,
    END_STAT        = END_STAT,
    NPKID           = NPKID,
    ABO_A           = ABO_A,
    ABO_B           = ABO_B,
    ABO_AB          = ABO_AB,
    ABO_O           = ABO_O,
    ETH_Asian       = ETH_Asian,
    ETH_Black       = ETH_Black,
    ETH_Hispanic    = ETH_Hispanic,
    ETH_White       = ETH_White,
    DWFG_KI         = DWFG_KI,
    GSTATUS_KI      = GSTATUS_KI,
    HIST_CIG_DON    = HIST_CIG_DON,
    SBP             = SBP,
    GENDER_DON      = GENDER_DON,
    B1              = B1,
    DB1             = DB1,
    B2              = B2,
    DB2             = DB2,
    DDR1            = DDR1,
    DR1             = DR1,
    DDR2            = DDR2,
    DR2             = DR2,
    HLA_DR_mismatches = HLA_DR_mismatches,
    HLA_B_mismatches  = HLA_B_mismatches,
    ABO_incompatible  = ABO_incompatible,
    CITIZENSHIP       = CITIZENSHIP,
    PRI_PAYMENT_TCR_KI = PRI_PAYMENT_TCR_KI,
    LISTING_CTR_CODE  = LISTING_CTR_CODE,
    WORK_INCOME_TRR   = WORK_INCOME_TRR,
    WORK_INCOME_TCR   = WORK_INCOME_TCR,
    TRTREJ1Y_KI       = TRTREJ1Y_KI,
    TRTREJ6M_KI       = TRTREJ6M_KI,
    FUNC_STAT_TRF     = FUNC_STAT_TRF
  ) %>%
  tidyr::drop_na()


dim(dmu_data)
table(dmu_data$Group)/nrow(dmu_data)

write.csv(dmu_data, "dmu_data.csv", row.names = F) # For Simulation in Supplementary Materials (S2) 

####################################################################################################
#########################################  RESAMPLING (Section 3.2) ################################
####################################################################################################

## ESRD Prevalence Counts by Year (Table B.1, USRDS 2024)
esrd_counts_by_year <- list(
  "2010" = c(Asian = 25048, Black = 183948, Hispanic = 94650, White = 276777),
  "2011" = c(Asian = 26676, Black = 190914, Hispanic = 100440, White = 282207),
  "2012" = c(Asian = 28345, Black = 197394, Hispanic = 105730, White = 289956),
  "2013" = c(Asian = 30148, Black = 204333, Hispanic = 111778, White = 298822),
  "2014" = c(Asian = 31956, Black = 211315, Hispanic = 117768, White = 308436),
  "2015" = c(Asian = 33800, Black = 217563, Hispanic = 124270, White = 318856),
  "2016" = c(Asian = 35825, Black = 222328, Hispanic = 130416, White = 328987),
  "2017" = c(Asian = 37810, Black = 226188, Hispanic = 136111, White = 338207),
  "2018" = c(Asian = 40088, Black = 230048, Hispanic = 142358, White = 348471),
  "2019" = c(Asian = 42642, Black = 235844, Hispanic = 149856, White = 356921)
)

## Convert counts to proportions
get_proportions <- function(counts) {
  total <- sum(counts)
  round(counts / total, 5)
}

## Stratified sampling within one year
stratified_sample <- function(data, group_var, desired_props, total_samples) {
  group_counts   <- table(data[[group_var]])
  desired_counts <- round(desired_props[names(group_counts)] * total_samples)
  
  split(data, data[[group_var]]) |>
    purrr::map2_dfr(
      .y = names(desired_counts),
      ~ dplyr::slice_sample(.x, n = min(nrow(.x), desired_counts[.y]))
    )
}

## Drop unused factor levels (and convert characters → factors)
drop_all_levels <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.character(col)) col <- factor(col)
    if (is.factor(col))   col <- droplevels(col)
    col
  })
  df
}

## Main resampling wrapper
resample_to_esrd_mix <- function(dmu_data, esrd_counts_by_year, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  years <- intersect(
    names(esrd_counts_by_year),
    as.character(unique(dmu_data$Year))
  )
  
  balanced_data_list <- lapply(years, function(year) {
    year_data <- dmu_data |> dplyr::filter(Year == as.integer(year))
    if (!nrow(year_data)) return(NULL)
    
    desired_props <- get_proportions(esrd_counts_by_year[[year]])
    group_sizes   <- table(year_data$Group)
    
    # Max feasible total sample while honoring desired proportions
    max_possible  <- group_sizes / desired_props[names(group_sizes)]
    total_samples <- floor(min(max_possible))
    
    stratified_sample(
      data          = year_data,
      group_var     = "Group",
      desired_props = desired_props,
      total_samples = total_samples
    )
  })
  
  dplyr::bind_rows(balanced_data_list)
}

DEA_data <- resample_to_esrd_mix(dmu_data, esrd_counts_by_year, seed = 14) %>%
  drop_all_levels() %>%
  group_by(Year) %>%
  mutate(
    Priority_Score = WaitlistDuration - mean(WaitlistDuration, na.rm = TRUE),
    Access_Score   = LKDPI            - mean(LKDPI, na.rm = TRUE),
    Outcome_Score  = GTIME_KI         - mean(GTIME_KI, na.rm = TRUE)
  ) %>%
  ungroup()



write.csv(DEA_data, "DEA_data.csv", row.names = F)


################################################################################
############################## Data Summaries ##################################
################################################################################

################################################################################
# Data Summaries for Fairness Criterion 
################################################################################

################################################################################
#                                    Raw Scores
################################################################################

# Raw Waitlist time Summaries
DEA_data %>%
  group_by(Group) %>%
  summarise(
    Mean = formatC(mean(WaitlistDuration, na.rm = TRUE), format = "f", digits = 1),
    SD = formatC(sd(WaitlistDuration, na.rm = TRUE), format = "f", digits = 1),
    Median = formatC(median(WaitlistDuration, na.rm = TRUE), format = "f", digits = 2),
    Q1 = formatC(quantile(WaitlistDuration, 0.25, na.rm = TRUE), format = "f", digits = 2),
    Q3 = formatC(quantile(WaitlistDuration, 0.75, na.rm = TRUE), format = "f", digits = 2)
  )

# Raw LKDPI Summaries
DEA_data %>%
  group_by(Group) %>%
  summarise(
    Mean = formatC(mean(LKDPI, na.rm = TRUE), format = "f", digits = 2),
    SD = formatC(sd(LKDPI, na.rm = TRUE), format = "f", digits = 1),
    Median = formatC(median(LKDPI, na.rm = TRUE), format = "f", digits = 2),
    Q1 = formatC(quantile(LKDPI, 0.25, na.rm = TRUE), format = "f", digits = 2),
    Q3 = formatC(quantile(LKDPI, 0.75, na.rm = TRUE), format = "f", digits = 2)
  )


# Raw Graft Lifespan Summaries
DEA_data %>%
  group_by(Group) %>%
  summarise(
    Mean = formatC(mean(GTIME_KI, na.rm = TRUE), format = "f", digits = 1),
    SD = formatC(sd(GTIME_KI, na.rm = TRUE), format = "f", digits = 1),
    Median = formatC(median(GTIME_KI, na.rm = TRUE), format = "f", digits = 2),
    Q1 = formatC(quantile(GTIME_KI, 0.25, na.rm = TRUE), format = "f", digits = 2),
    Q3 = formatC(quantile(GTIME_KI, 0.75, na.rm = TRUE), format = "f", digits = 2)
  )


################################################################################
#                               Relative Scores
################################################################################

# Relative Waitlist time Summaries
DEA_data %>%
  group_by(Group) %>%
  summarise(
    Mean = formatC(mean(Priority_Score, na.rm = TRUE), format = "f", digits = 1),
    SD = formatC(sd(Priority_Score, na.rm = TRUE), format = "f", digits = 1),
    Median = formatC(median(Priority_Score, na.rm = TRUE), format = "f", digits = 2),
    Q1 = formatC(quantile(Priority_Score, 0.25, na.rm = TRUE), format = "f", digits = 2),
    Q3 = formatC(quantile(Priority_Score, 0.75, na.rm = TRUE), format = "f", digits = 2)
  )

# Relative LKDPI Summaries
DEA_data %>%
  group_by(Group) %>%
  summarise(
    Mean = formatC(mean(Access_Score, na.rm = TRUE), format = "f", digits = 1),
    SD = formatC(sd(Access_Score, na.rm = TRUE), format = "f", digits = 1),
    Median = formatC(median(Access_Score, na.rm = TRUE), format = "f", digits = 2),
    Q1 = formatC(quantile(Access_Score, 0.25, na.rm = TRUE), format = "f", digits = 2),
    Q3 = formatC(quantile(Access_Score, 0.75, na.rm = TRUE), format = "f", digits = 2)
  )


# Relative Graft Lifespan Summaries
DEA_data %>%
  group_by(Group) %>%
  summarise(
    Mean = formatC(mean(Outcome_Score, na.rm = TRUE), format = "f", digits = 1),
    SD = formatC(sd(Outcome_Score, na.rm = TRUE), format = "f", digits = 1),
    Median = formatC(median(Outcome_Score, na.rm = TRUE), format = "f", digits = 2),
    Q1 = formatC(quantile(Outcome_Score, 0.25, na.rm = TRUE), format = "f", digits = 2),
    Q3 = formatC(quantile(Outcome_Score, 0.75, na.rm = TRUE), format = "f", digits = 2)
  )


################################################################################
# Additional summaries:
DEA_data %>%
  group_by(Year, Group) %>%
  summarise(Average_Score = mean(Priority_Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(
    names_from = Year,
    values_from = Average_Score
  )


DEA_data %>%
  group_by(Year, Group) %>%
  summarise(Average_Score = mean(Access_Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(
    names_from = Year,
    values_from = Average_Score
  )

DEA_data %>%
  group_by(Year, Group) %>%
  summarise(Average_Score = mean(Outcome_Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(
    names_from = Year,
    values_from = Average_Score
  )

################################################################################
################################################################################


################################################################################
######################### Data Imbalance Plot (Fig. 2) #########################
################################################################################

plot_year <- 2019

# Ensure consistent group ordering
DEA_data <- DEA_data %>%
  dplyr::mutate(
    Group = factor(Group, levels = c("Asian", "Black", "Hispanic", "White"))
  )

dmu_data <- dmu_data %>%
  dplyr::mutate(
    Group = factor(Group, levels = c("Asian", "Black", "Hispanic", "White"))
  )

# ----- Original (imbalanced) distribution in the transplant data -----
original_data <- dmu_data %>%
  dplyr::filter(Year == plot_year) %>%
  dplyr::count(Group, name = "n") %>%
  dplyr::mutate(
    Percentage = n / sum(n),
    Dataset    = "Imbalanced Data"
  )

# ----- Target ESRD prevalence distribution (from registry) -----
esrd_props_2019 <- get_proportions(esrd_counts_by_year[["2019"]])

resampled_data <- DEA_data %>%
  dplyr::filter(Year == plot_year) %>%
  dplyr::count(Group, name = "n") %>%
  dplyr::mutate(
    Percentage = esrd_props_2019[as.character(Group)],
    Dataset    = "ESRD Prevalence"
  )

# ----- Combine & factorise dataset label -----
combined_data <- dplyr::bind_rows(original_data, resampled_data) %>%
  dplyr::mutate(
    Dataset = factor(Dataset, levels = c("Imbalanced Data", "ESRD Prevalence"))
  )

# ----- Colors and plot -----
color_palette <- c(
  "Asian"    = "#08306b",
  "Black"    = "#1f78b4",
  "Hispanic" = "#6baed6",
  "White"    = "#b3cde3"
)

ggplot(combined_data, aes(x = Dataset, y = Percentage, fill = Group)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(width = 0.8),
    width    = 0.8
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", Percentage * 100)),
    position = position_dodge(width = 0.8),
    vjust    = -0.5,
    size     = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = color_palette) +
  labs(x = "", y = "Relative Frequency") +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1),
    limits = c(0, 0.75)
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold"),
    axis.text.y  = element_text(size = 16),
    axis.title.y = element_text(size = 20, face = "bold",
                                margin = ggplot2::margin(r = 10)),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16)
  )
