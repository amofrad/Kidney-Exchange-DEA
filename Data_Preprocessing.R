library(readr)
library(tidyverse)
library(dplyr)
library(lubridate)
library(ggplot2)
library(purrr)
library(forcats)


data <- read_delim("KIDPAN_DATA.DAT", 
                   delim = "\t", escape_double = FALSE,
                   col_names = FALSE, trim_ws = TRUE)

colnames(data)[1:265] <- c("WL_ORG", "COD_WL", "COD_OSTXT_WL", "NUM_PREV_TX", "CURRENT_PRA", "PEAK_PRA", "USE_WHICH_PRA", "CREAT_CLEAR", "GFR", "DONATION", "ON_DIALYSIS", "MAX_KDPI_LOCAL_ZERO_ABDR", "MAX_KDPI_LOCAL_NON_ZERO_ABDR", "MAX_KDPI_IMPORT_ZERO_ABDR", "MAX_KDPI_IMPORT_NON_ZERO_ABDR", "C_PEPTIDE", "C_PEPTIDEDATE", "A2A2B_ELIGIBILITY", "A1", "A2", "B1", "B2", "DR1", "DR2", "ANTIBODY_TESTED", "GENDER", "ABO", "WGT_KG_TCR", "HGT_CM_TCR", "BMI_TCR", "CITIZENSHIP", "CITIZEN_COUNTRY", "PERM_STATE", "EDUCATION", "FUNC_STAT_TCR", "DGN_TCR", "DGN_OSTXT_TCR", "DGN2_TCR", "DGN2_OSTXT_TCR", "DIAB", "DRUGTRT_COPD", "TOT_SERUM_ALBUM", "C_PEPTIDE_PA_TCR", "HBA1C_PA_TCR", "SSDMF_DEATH_DATE", "INIT_CURRENT_PRA", "INIT_PEAK_PRA", "INIT_STAT", "INIT_WGT_KG", "INIT_HGT_CM", "INIT_CPRA", "END_CPRA", "INIT_EPTS", "END_EPTS", "REM_CD", "DAYSWAIT_CHRON", "END_STAT", "INIT_AGE","ACTIVATE_DATE", "CREAT_CLEAR_DATE", "DEATH_DATE", "DIALYSIS_DATE", "END_DATE", "GFR_DATE", "INIT_DATE", "WT_QUAL_DATE", "ETHNICITY", "ETHCAT", "PT_CODE","INIT_BMI_CALC", "END_BMI_CALC", "DAYSWAIT_ALLOC", "COMPOSITE_DEATH_DATE", "WLHR", "WLHL", "WLIN", "WLKI", "WLKP", "WLLI", "WLLU", "WLPA", "WLPI", "WLVC", "REGION", "INACT_REASON_CD", "BW4", "BW6", "C1", "C2", "DR51", "DR51_2", "DR52", "DR52_2", "DR53", "DR53_2", "DQ1", "DQ2", "WL_ID_CODE", "PERIP_VASC", "EXH_PERIT_ACCESS", "AGE_DIAB", "EXH_VASC_ACCESS", "YR_ENTRY_US_TCR", "WORK_INCOME_TCR", "ACADEMIC_PRG_TCR", "ACADEMIC_LEVEL_TCR", "MALIG_TCR_KI", "PRI_PAYMENT_TCR_KI", "MALIG_TCR_PA", "PRI_PAYMENT_TCR_PA", "PREV_TX", "PREV_KI_TX", "PREV_PA_TX", "ACADEMIC_LEVEL_TRR", "ACADEMIC_PRG_TRR", "FUNC_STAT_TRR", "MALIG_TRR", "MALIG_OSTXT_TRR", "MALIG_TY_TRR", "PERM_STATE_TRR", "PRI_PAYMENT_TRR_KI", "PRI_PAYMENT_CTRY_TRR_KI", "WORK_INCOME_TRR", "TX_DATE", "ACUTE_REJ_EPI_KI", "CREAT_TRR", "FIRST_WK_DIAL", "ORG_REC_ON", "PREV_PREG", "REC_ON_ICE", "REC_ON_PUMP", "SERUM_CREAT", "L_FIN_FLOW_RATE_TX", "L_FIN_RESIST_TX", "R_FIN_FLOW_RATE_TX", "R_FIN_RESIST_TX", "PRE_TX_TXFUS", "FIN_RESIST_TX", "TXHRT", "TXINT", "TXKID", "TXLIV", "TXLNG", "TXPAN", "TXVCA", "PRI_PAYMENT_CTRY_TCR_KI", "PREV_MALIG_TY", "PREV_MALIG_TY_OSTXT", "RDA1", "RDA2", "RDB1", "RDB2", "RDDR1", "RDDR2", "DON_RETYP", "RESUM_MAINT_DIAL_DT", "DA1", "DA2", "DB1", "DB2", "DDR1", "DDR2", "RA1", "RA2", "RB1", "RB2", "RDR1", "RDR2", "AMIS", "BMIS", "DRMIS", "HLAMIS", "NPKID", "NPPAN", "END_CPRA_DETAIL", "HAPLO_TY_MATCH_DON", "AGE_DON", "DDAVP_DON", "CMV_OLD_LIV_DON", "CMV_DON", "CMV_TEST_DON", "EBV_TEST_DON", "HBV_TEST_DON", "HCV_TEST_DON", "CMV_NUCLEIC_DON", "CMV_IGG_DON", "CMV_IGM_DON", "EBV_DNA_DON", "EBV_IGG_DON", "EBV_IGM_DON", "HBV_CORE_DON", "HBV_SUR_ANTIGEN_DON", "ETHCAT_DON", "COD_CAD_DON", "DEATH_CIRCUM_DON", "DEATH_MECH_DON", "CITIZENSHIP_DON", "HEP_C_ANTI_DON", "HCV_RNA_DON", "ABO_DON", "DON_TY", "GENDER_DON", "HOME_STATE_DON", "WARM_ISCH_TM_DON", "HCV_RIBA_DON", "HCV_ANTIBODY_DON", "LIV_DON_TY", "CITIZEN_COUNTRY_DON", "COD_OSTXT_DON", "CONTROLLED_DON", "CORE_COOL_DON", "NON_HRT_DON", "ANTIHYPE_DON", "BLOOD_INF_DON", "BLOOD_INF_CONF_DON", "BUN_DON", "CREAT_DON", "DOBUT_DON_OLD", "DOPAMINE_DON_OLD", "HTLV1_OLD_DON", "HTLV2_OLD_DON", "OTH_DON_MED1_OSTXT_DON_OLD", "OTH_DON_MED2_OSTXT_DON_OLD", "OTH_DON_MED3_OSTXT_DON_OLD", "OTHER_INF_DON", "OTHER_INF_CONF_DON", "OTHER_INF_OSTXT_DON", "PRETREAT_MED_DON_OLD", "PT_DIURETICS_DON", "PT_STEROIDS_DON", "PT_T3_DON", "PT_T4_DON", "PT_OTH2_OSTXT_DON", "PT_OTH3_OSTXT_DON", "PT_OTH4_OSTXT_DON", "PT_OTH1_OSTXT_DON", "PULM_INF_DON", "PULM_INF_CONF_DON", "SGOT_DON", "SGPT_DON", "TBILI_DON", "URINE_INF_DON", "URINE_INF_CONF_DON", "VASODIL_DON", "VDRL_DON", "CLIN_INFECT_DON", "HYPERTENS_DUR_DON", "CANCER_FREE_INT_DON", "CANCER_OTH_OSTXT_DON", "CONTIN_ALCOHOL_OLD_DON", "CONTIN_CIG_DON", "CONTIN_IV_DRUG_OLD_DON", "CONTIN_COCAINE_DON", "CONTIN_OTH_DRUG_DON", "DIET_DON", "DIURETICS_DON", "EXTRACRANIAL_CANCER_DON"         , "HIST_ALCOHOL_OLD_DON", "CANCER_SITE_DON", "HIST_CIG_DON", "DIABDUR_DON", "HIST_COCAINE_DON", "HIST_HYPERTENS_DON", "HIST_IV_DRUG_OLD_DON", "INSULIN_DEP_DON")
colnames(data)[266:491] <- c("INTRACRANIAL_CANCER_DON", "OTHER_HYPERTENS_MED_DON", "HIST_CANCER_DON", "HIST_INSULIN_DEP_DON", "INSULIN_DUR_DON", "HIST_DIABETES_DON", "HIST_OTH_DRUG_DON", "SKIN_CANCER_DON", "DIABETES_DON", "LIV_DON_TY_OSTXT", "HEPARIN_DON", "ARGININE_DON", "INSULIN_DON", "HGT_CM_DON_CALC", "WGT_KG_DON_CALC", "BMI_DON_CALC", "KDPI", "KDRI_MED", "KDRI_RAO", "HBV_NAT_DON", "HCV_NAT_DON", "HIV_NAT_DON", "END_STAT_KI", "CREAT6M", "CREAT1Y", "DIAL_DATE", "RETXDATE_KI", "FAILDATE_KI", "PUMP_KI", "ABO_MAT", "AGE", "DISTANCE", "RESUM_MAINT_DIAL", "DIAL_TRR", "DIAG_KI", "DIAG_OSTXT_KI", "COLD_ISCH_KI", "GRF_STAT_KI", "GRF_FAIL_CAUSE_OSTXT_KI", "GRF_FAIL_CAUSE_TY_KI", "DWFG_KI", "PRVTXDIF_KI", "GTIME_KI", "GSTATUS_KI", "COD_KI", "COD_OSTXT_KI", "COD2_KI", "COD2_OSTXT_KI", "COD3_KI", "COD3_OSTXT_KI", "DAYSWAIT_CHRON_KI", "TX_PROCEDUR_TY_KI", "TRTREJ1Y_KI", "TRTREJ6M_KI", "MULTIORG", "PRI_PAYMENT_TRR_PA", "PRI_PAYMENT_CTRY_TRR_PA", "ART_RECON", "ART_RECON_OSTXT", "DUCT_MGMT", "DUCT_MGMT_OSTXT", "GRF_PLACEM", "PRE_AVG_INSULIN_USED_TRR", "PRE_AVG_INSULIN_USED_OLD_TRR", "ACUTE_REJ_EPI_PA", "PA_PRESERV_TM", "VASC_MGMT", "VEN_EXT_GRF", "INSULIN_PA", "INSULIN_RESUMED_DATE_PA", "INSULIN_DOSAGE_PA", "INSULIN_DURATION_PA", "METHOD_BLOOD_SUGAR_CONTROL_PA", "BLOOD_SUGAR_MEDICATION_PA", "BLOOD_SUGAR_MED_RESUMED_DATE_PA", "BLOOD_SUGAR_DIET_PA", "C_PEPTIDE_PA_TRR", "HBA1C_PA_TRR", "INSULIN_DOSAGE_OLD_PA", "PRI_PAYMENT_CTRY_TCR_PA", "PK_DA1", "PK_DA2", "PK_DB1", "PK_DB2", "PK_DDR1", "PK_DDR2", "ENTERIC_DRAIN", "ENTERIC_DRAIN_DT", "END_STAT_PA", "FAILDATE_PA", "DIAG_PA", "DIAG_OSTXT_PA", "GRF_STAT_PA", "GRF_FAIL_CAUSE_OSTXT_PA", "GRF_FAIL_CAUSE_TY_PA", "OTH_GRF_FAIL_CAUSE_OSTXT_PA", "GRF_VASC_THROMB_PA", "INFECT_PA", "BLEED_PA", "ANAST_LK_PA", "REJ_ACUTE_PA", "REJ_HYPER_PA", "BIOP_ISLET_PA", "PANCREATIT_PA", "REJ_CHRONIC_PA", "PX_NON_COMPL_PA", "RETXDATE_PA", "PRVTXDIF_PA", "GTIME_PA", "GSTATUS_PA", "COD_PA", "COD_OSTXT_PA", "COD2_PA", "COD2_OSTXT_PA", "COD3_PA", "COD3_OSTXT_PA", "DAYSWAIT_CHRON_PA", "TX_PROCEDUR_TY_PA", "TRTREJ1Y_PA", "TRTREJ6M_PA", "ORGAN", "CMV_IGG", "CMV_IGM", "EBV_SEROSTATUS", "HBV_CORE", "HBV_SUR_ANTIGEN", "HCV_SEROSTATUS", "HIV_SEROSTATUS", "CMV_STATUS", "HBV_SURF_TOTAL", "HIV_NAT", "HCV_NAT", "HBV_NAT", "PREV_TX_ANY", "PREV_TX_ANY_N", "TX_TYPE", "MED_COND_TRR", "PX_STAT", "PX_STAT_DATE", "PREV_KI_DATE", "FUNC_STAT_TRF", "SHARE_TY", "PSTATUS", "PTIME", "LOS", "PAYBACK", "ECD_DONOR", "AGE_GROUP", "MALIG", "MALIG_TY_OSTXT", "MALIG_TY", "HGT_CM_CALC", "WGT_KG_CALC", "BMI_CALC", "STATUS_TCR", "STATUS_TRR", "STATUS_DDR", "VAL_DT_DDR", "STATUS_LDR", "VAL_DT_LDR", "VAL_DT_TCR", "VAL_DT_TRR", "LT_ONE_WEEK_DON", "REJ_BIOPSY", "REJCNF_KI", "REJTRT_KI", "REJCNF_PA", "REJTRT_PA", "TRR_ID_CODE", "ADMISSION_DATE", "DISCHARGE_DATE", "COMPL_ABSC", "COMPL_ANASLK", "COMPL_PANCREA", "OTH_COMPL_OSTXT", "SURG_INCIS", "OPER_TECH", "EDUCATION_DON", "KI_CREAT_PREOP", "KI_PROC_TY", "PRI_PAYMENT_DON", "PRI_PAYMENT_CTRY_DON", "MEDICARE_DON", "MEDICAID_DON", "OTH_GOVT_DON", "PRIV_INS_DON", "HMO_PPO_DON", "SELF_DON", "DONATION_DON", "FREE_DON", "RECOV_OUT_US", "RECOV_COUNTRY", "PROTEIN_URINE", "LIPASE", "AMYLASE", "INOTROP_AGENTS", "CARDARREST_NEURO", "RESUSCIT_DUR", "INOTROP_SUPPORT_DON", "TATTOOS", "LT_KI_BIOPSY", "LT_KI_GLOMERUL", "RT_KI_BIOPSY", "RT_KI_GLOMERUL", "REFERRAL_DATE", "RECOVERY_DATE", "ADMIT_DATE_DON", "DONOR_ID", "HBSAB_DON", "EBV_IGG_CAD_DON", "EBV_IGM_CAD_DON", "HBV_DNA_DON", "CDC_RISK_HIV_DON", "INO_PROCURE_AGENT_1", "INO_PROCURE_AGENT_2", "INO_PROCURE_AGENT_3", "INO_PROCURE_OSTXT_1", "INO_PROCURE_OSTXT_2", "INO_PROCURE_OSTXT_3", "DATA_TRANSPLANT", "DATA_WAITLIST", "CTR_CODE", "OPO_CTR_CODE", "INIT_OPO_CTR_CODE", "END_OPO_CTR_CODE", "LISTING_CTR_CODE")

# KPD only
data <- data %>% filter((LIV_DON_TY == 9))

write.csv(data, "KID_LIVING.csv", col.names = F)

################################################################################
################################# DATA PREP ####################################
################################################################################

kidney_data <- read.csv("KID_LIVING.csv")

kidney_data$TX_YEAR <- year(mdy(kidney_data$TX_DATE))
kidney_data <- kidney_data %>% filter(TX_YEAR >= 2010, TX_YEAR <= 2019)

kidney_data$DAYSWAIT_CHRON_KI <- ifelse(kidney_data$DAYSWAIT_CHRON_KI == ".", 
                                        NA, 
                                        kidney_data$DAYSWAIT_CHRON_KI)

kidney_data$DAYSWAIT_CHRON_KI <- as.numeric(kidney_data$DAYSWAIT_CHRON_KI)

kidney_data$DIAL_TIME <- as.numeric(
  ifelse(kidney_data$DIALYSIS_DATE == ".",
         0,
         as.Date(kidney_data$INIT_DATE, "%m/%d/%Y") - as.Date(kidney_data$DIALYSIS_DATE, "%m/%d/%Y"))
)

kidney_data <- kidney_data %>% 
  mutate(ETHCAT = case_when(
    ETHCAT == 1 ~ "White",
    ETHCAT == 2 ~ "Black",
    ETHCAT == 4 ~ "Hispanic",
    ETHCAT == 5 ~ "Asian",
    ETHCAT == 6 ~ "Other",
    ETHCAT == 7 ~ "Other",
    ETHCAT == 9 ~ "Other",
    ETHCAT == 998 ~ "Other",
    TRUE ~ as.character(ETHCAT) 
  ),
  ETHCAT_DON = case_when(
    ETHCAT_DON == 1 ~ "White",
    ETHCAT_DON == 2 ~ "Black",
    ETHCAT_DON == 4 ~ "Hispanic",
    ETHCAT_DON == 5 ~ "Asian",
    ETHCAT_DON == 6 ~ "Other",
    ETHCAT_DON == 7 ~ "Other",
    ETHCAT_DON == 9 ~ "Other",
    ETHCAT_DON == 998 ~ "Other",
    TRUE ~ as.character(ETHCAT_DON) 
  )
  )

kidney_data$ETH_White <- as.numeric(kidney_data$ETHCAT == "White")
kidney_data$ETH_Asian <- as.numeric(kidney_data$ETHCAT == "Asian")
kidney_data$ETH_Black <- as.numeric(kidney_data$ETHCAT == "Black")
kidney_data$ETH_Hispanic <- as.numeric(kidney_data$ETHCAT == "Hispanic")
kidney_data$ETHCAT <- as.factor(kidney_data$ETHCAT)
kidney_data$ETHCAT_DON <- as.factor(kidney_data$ETHCAT_DON)

Donors <- read_delim("LIVING_DONOR_DATA.DAT", 
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE)

colnames(Donors) <- c("REOP_BILIARY", "REOP_BILIARY_DT", "DON_DATE", "AGE_DON", "ETHCAT_DON", "REGION", "LUNG_RECOV", "KIDNEY_RECOV", "LIVER_RECOV", "DON_ORG2", "CITIZENSHIP", "LIV_DON_TY", "LIV_DON_TY_OSTXT", "EDUCATION", "KI_CREAT_PREOP", "BP_PREOP_SYST", "BP_PREOP_DIAST", "COD", "KI_PROC_TY", "MARITAL_STAT", "HEALTH_INS", "FUNC_STAT", "PHYSICAL_CAPACITY", "WORK_INCOME", "HIST_CANCER", "HIST_CANCER_OSTXT", "CANCER_FREE", "COD_OSTXT", "HIST_HYPER", "HYPER_DIET", "HYPER_DIUR", "HYPER_MEDS", "PREOP_URINE_RATIO", "DIABETES", "MACRO_FAT", "MICRO_FAT", "CONVERT_OPEN_KI", "NON_AUTO_BLOOD", "PRBC_UNITS", "PLATELETS_UNITS", "FFP_UNITS", "VASC_COMP_KI", "VASC_COMP_KI_INTER", "VASC_COMP_KI_INTER_OSTXT", "OTH_COMP_KI", "OTH_COMP_KI_INTER", "OTH_COMP_KI_INTER_OSTXT", "REOPERATION_KI", "REOP_BLEED_KI", "REOP_HERNIA_KI", "REOP_BOWEL_KI", "REOP_VASC_KI", "REOP_OTH_KI", "REOP_OTH_KI_OSTXT", "READMISSION_KI", "READMISSION_KI_REASON", "READMISSION_KI_OSTXT", "OTH_INTER_PROC_KI", "OTH_INTER_PROC_KI_OSTXT", "KI_CREAT_POSTOP", "BP_POSTOP_SYST", "BP_POSTOP_DIAST", "HYPERTENSION", "POSTOP_URINE_RATIO", "PREOP_BILI", "PREOP_SGOT_AST", "PREOP_ALK_PHOS", "PREOP_ALBUM", "PREOP_INR", "BIOPSY_LI", "LI_PROC_TY", "BILIARY_COMP", "BILIARY_COMP_GRADE", "VASC_COMP_LI", "VASC_COMP_LI_INTER", "VASC_COMP_LI_INTER_OSTXT", "OTH_COMP_LI", "OTH_COMP_LI_INTER", "OTH_COMP_LI_INTER_OSTXT", "REOPERATION_LI", "REOP_BLEED_LI", "REOP_HERNIA_LI", "REOP_BOWEL_LI", "REOP_VASC_LI", "REOP_OTH_LI", "REOP_OTH_LI_OSTXT", "REOP_LI_FAIL", "READMISSION_LI", "READMISSION_LI_REASON", "READMISSION_LI_OSTXT", "OTH_INTER_PROC_LI", "OTH_INTER_PROC_LI_OSTXT", "POSTOP_SGOT_AST", "POSTOP_ALK_PHOS", "POSTOP_ALBUM", "POSTOP_CREAT_LI", "POSTOP_INR", "POSTOP_SGPT_ALT", "POSTOP_BILI", "PREOP_FVC_BEFORE", "PREOP_FVC_AFTER", "PREOP_FEV1_BEFORE", "PREOP_FEV1_AFTER", "PREOP_FEF_BEFORE", "PREOP_FEF_AFTER", "PREOP_TLC_BEFORE", "PREOP_TLC_AFTER", "PREOP_LUNG_CAP", "PREOP_PAO2", "HIST_CIG", "PACK_YRS", "DUR_ABSTINENCE", "TOBACCO_USE", "LU_PROC_TY", "CONVERT_OPEN_LU", "INTRAOP_COMP", "INTRAOP_COMP_REASON", "SACRIFICE_LOBE", "ARRHYTHMIA", "ANESTHETIC_COMP", "INTRAOP_COMP_OSTXT", "LU_COMP", "LU_COMP_REASON", "LU_COMP_OSTXT", "THORAC_TUBES", "ARRHYTHMIA_POSTOP", "READMISSION_LU", "READMISSION_LU_REASON", "READMISSION_LU_OSTXT", "PREDON_HGT", "PREDON_WGT", "PREOP_SGPT_ALT", "PREOP_CREAT_LI", "PREOP_URINE_PROTEIN", "POSTOP_URINE_PROTEIN", "PX_STAT", "CITIZEN_COUNTRY", "DEATH_DT", "ORG_RECOVERY_DT", "INIT_DISCHARGE_DT", "REOP_BLEED_KI_DT", "REOP_HERNIA_KI_DT", "REOP_BOWEL_KI_DT", "REOP_VASC_KI_DT", "REOP_OTH_KI_DT", "READMISSION_KI_DT", "OTH_INTER_PROC_KI_DT", "POSTOP_TEST_DT", "REOP_BLEED_LI_DT", "REOP_HERNIA_LI_DT", "REOP_BOWEL_LI_DT", "REOP_VASC_LI_DT", "REOP_OTH_LI_DT", "REOP_LI_FAIL_DT", "READMISSION_LI_DT", "OTH_INTER_PROC_LI_DT", "READMISSION_LU_DT", "CMV_IGG", "CMV_IGM", "CMV_NUCLEIC", "EBV_IGG", "EBV_IGM", "HBV_CORE", "HBV_DNA", "HBV_SUR_ANTIGEN", "HCV_ANTIBODY", "HCV_RIBA", "HCV_RNA", "VIRUSES_TESTED", "CMV_TOTAL", "EBV_TOTAL", "HOME_STATE", "GENDER", "ABO", "YR_ENTRY_US", "WGT_KG", "DON_ORG", "STATUS_LDR", "VAL_DT_LDR", "DBW4", "DBW6", "DC1", "DC2", "DDP1", "DDP2", "DDPA1", "DDPA2", "DDR51", "DDR51_2", "DDR52", "DDR52_2", "DDR53", "DDR53_2", "DDQ1", "DDQ2", "DDQA1", "DDQA2", "RECOV_FACILITY_CODE", "DONOR_ID")

Donors$KI_CREAT_PREOP <- as.numeric(Donors$KI_CREAT_PREOP)
Donors$African_American <- Donors$ETHCAT_DON == 2

# Function to calculate eGFR
calculate_eGFR <- function(KI_CREAT_PREOP, AGE_DON, GENDER_DON, African_American) {
  kappa <- ifelse(GENDER_DON == "F", 0.7, 0.9)
  alpha <- ifelse(GENDER_DON == "F", -0.329, -0.411)
  eGFR <- 141 * min(KI_CREAT_PREOP/kappa, 1)^alpha * max(KI_CREAT_PREOP/kappa, 1)^-1.209 * 0.993^AGE_DON
  
  if(GENDER_DON == "F") {
    eGFR <- eGFR * 1.018
  }
  
  if(African_American) {
    eGFR <- eGFR * 1.159
  }
  
  return(eGFR)
}

# Keeping all rows from kidney_data (left join)
merged_data <- merge(kidney_data, Donors, by = "DONOR_ID", all.x = TRUE)

# Recreate variables after merging
age <- as.numeric(merged_data$AGE_DON.x)
eGFR <- numeric(nrow(merged_data))
for (i in 1:nrow(merged_data)) {
  if (complete.cases(merged_data$KI_CREAT_PREOP.y[i], merged_data$AGE_DON.y[i], merged_data$GENDER.y[i], merged_data$African_American[i])) {
    eGFR[i] <- calculate_eGFR(merged_data$KI_CREAT_PREOP.y[i], merged_data$AGE_DON.y[i], merged_data$GENDER.y[i], merged_data$African_American[i])
  } else {
    eGFR[i] <- NA
  }
}

merged_data$eGFR <- eGFR

BMI <- as.numeric(merged_data$BMI_DON_CALC)
african_american <- merged_data$ETHCAT_DON.x == "Black" 
history_of_cigarette_use <- merged_data$HIST_CIG_DON == "Y" 
SBP <- as.numeric(merged_data$BP_PREOP_SYST)
both_male <- (merged_data$GENDER.x == "M") & (merged_data$GENDER_DON == "M")
ABO_incompatible <- merged_data$ABO_MAT == 3
unrelated <- !(as.factor(merged_data$LIV_DON_TY.x) %in% 1:6)
HLA_B_mismatches <- ((merged_data$B1) == (merged_data$DB1)) + ((merged_data$B2) == (merged_data$DB2))
HLA_DR_mismatches <- ((merged_data$DDR1) == (merged_data$DR1)) + ((merged_data$DDR2) == (merged_data$DR2))
D_RWR <- as.numeric(merged_data$WGT_KG_DON_CALC) / as.numeric(merged_data$WGT_KG_CALC)



calculate_LKDPI <- function(age, eGFR, BMI, african_american, history_of_cigarette_use, SBP, both_male, ABO_incompatible, unrelated, HLA_B_mismatches, HLA_DR_mismatches, D_RWR) {
  # Initialize the LKDPI with the constant term
  LKDPI <- -11.30
  
  # Add age component if age is over 50
  if(age > 50) {
    LKDPI <- LKDPI + 1.85 * (age - 50)
  }
  
  LKDPI <- LKDPI - 0.381 * eGFR + 1.17 * BMI
  
  if(african_american) {
    LKDPI <- LKDPI + 22.34
  }
  
  if(history_of_cigarette_use) {
    LKDPI <- LKDPI + 14.33
  }
  
  LKDPI <- LKDPI + 0.44 * SBP
  
  if(both_male) {
    LKDPI <- LKDPI - 21.68
  }
  
  if(ABO_incompatible) {
    LKDPI <- LKDPI + 27.30
  }
  
  if(unrelated) {
    LKDPI <- LKDPI - 10.61
  }
  
  LKDPI <- LKDPI + 8.57 * HLA_B_mismatches + 8.26 * HLA_DR_mismatches
  
  # Adjust D/RWR component
  LKDPI <- LKDPI - 50.87 * min(D_RWR, 0.9)
  
  return(LKDPI)
}

###################

LKDPI <- numeric(nrow(merged_data))

for (i in 1:nrow(merged_data)) {
  if (complete.cases(age[i], eGFR[i], BMI[i], african_american[i], history_of_cigarette_use[i],
                     SBP[i], both_male[i], ABO_incompatible[i], unrelated[i],
                     HLA_B_mismatches[i], HLA_DR_mismatches[i], D_RWR[i])) {
    LKDPI[i] <- calculate_LKDPI(age[i], eGFR[i], BMI[i], african_american[i],
                                history_of_cigarette_use[i], SBP[i], both_male[i],
                                ABO_incompatible[i], unrelated[i],
                                HLA_B_mismatches[i], HLA_DR_mismatches[i], D_RWR[i])
  } else {
    LKDPI[i] <- NA
  }
}

# Add LKDPI to the merged data
merged_data$LKDPI <- LKDPI
merged_data$SBP <- SBP
merged_data$HLA_B_mismatches <- HLA_B_mismatches
merged_data$HLA_DR_mismatches <- HLA_DR_mismatches
merged_data$ABO_incompatible <- ABO_incompatible

# Filter data to only keep observations with non-NA LKDPI
final_data <- merged_data[!is.na(merged_data$LKDPI), ]

names(final_data) <- sub("\\.x$", "", names(final_data))


vars_to_keep <- names(kidney_data)
vars_to_keep <- c(vars_to_keep, "LKDPI", "eGFR", "SBP", "HLA_B_mismatches", "HLA_DR_mismatches", "ABO_incompatible")
final_data <- final_data[, vars_to_keep]
kidney_data <- final_data

kidney_data <- kidney_data %>%
  filter(!is.na(ON_DIALYSIS)) %>%
  mutate(ON_DIALYSIS = ifelse(ON_DIALYSIS == "Y", 1, 0))
kidney_data$GTIME_KI <- as.numeric(kidney_data$GTIME_KI)
kidney_data$ON_DIALYSIS <- as.factor(kidney_data$ON_DIALYSIS)
kidney_data$ABO <- as.factor(kidney_data$ABO)
kidney_data$AGE <- as.numeric(kidney_data$AGE)
kidney_data$GENDER <- as.factor(kidney_data$GENDER)
kidney_data$PRA <- as.numeric(kidney_data$END_CPRA)
kidney_data$WGT_KG_CALC <- as.numeric(kidney_data$WGT_KG_CALC)
kidney_data$BMI <- as.numeric(kidney_data$BMI_CALC)
kidney_data$DIABETES_DON <- as.factor(kidney_data$DIABETES_DON)

kidney_data <- kidney_data %>% filter(ETHCAT %in% c("Asian", "Black", "Hispanic", "White"))
kidney_data$ETHCAT <- droplevels(kidney_data$ETHCAT)
kidney_data$ETHCAT_DON <- droplevels((kidney_data$ETHCAT_DON))

kidney_data <- kidney_data %>%
  mutate(
    PRI_PAYMENT_TRR_KI = case_when(
      PRI_PAYMENT_TRR_KI %in% 1:14 ~ as.factor(PRI_PAYMENT_TRR_KI),
      TRUE ~ as.factor("Other")  # Instead of NA
    ),
    ACUTE_REJ_EPI_KI = case_when(
      ACUTE_REJ_EPI_KI %in% 1:3 ~ as.factor(ACUTE_REJ_EPI_KI),
      TRUE ~ as.factor("Unknown")  # Instead of NA
    ),
    EDUCATION = case_when(
      EDUCATION %in% c(1:6, 996, 998) ~ as.factor(EDUCATION),
      TRUE ~ as.factor("Unknown")
    )
  )

kidney_data$PRI_PAYMENT_TRR_KI <- droplevels(kidney_data$PRI_PAYMENT_TRR_KI)
kidney_data$ACUTE_REJ_EPI_KI <- droplevels(kidney_data$ACUTE_REJ_EPI_KI)
kidney_data$EDUCATION <- droplevels(kidney_data$EDUCATION)
kidney_data$AGE_DON = as.numeric(kidney_data$AGE_DON)
kidney_data$HGT_CM_DON_CALC = as.numeric(kidney_data$HGT_CM_DON_CALC)
kidney_data$WGT_KG_DON_CALC = as.numeric(kidney_data$WGT_KG_DON_CALC)
kidney_data$HIST_HYPERTENS_DON = as.factor(kidney_data$HIST_HYPERTENS_DON)
kidney_data$DISTANCE = as.numeric(kidney_data$DISTANCE)
kidney_data$PREV_KI_TX = as.factor(kidney_data$PREV_KI_TX)
kidney_data$HIST_CANCER_DON = as.factor(kidney_data$HIST_CANCER_DON)
kidney_data$HIST_CIG_DON = as.factor(kidney_data$HIST_CIG_DON)
kidney_data$DIAG_KI = as.factor(kidney_data$DIAG_KI)
kidney_data$ACADEMIC_LEVEL_TCR = as.factor(kidney_data$ACADEMIC_LEVEL_TCR)
kidney_data$DIAB = as.factor(kidney_data$DIAB)
kidney_data$EDUCATION = as.factor(kidney_data$EDUCATION)
kidney_data$ABO_MAT = as.factor(kidney_data$ABO_MAT)
kidney_data$PREV_TX_ANY = as.factor(kidney_data$PREV_TX_ANY)
kidney_data$REGION = as.factor(kidney_data$REGION)
kidney_data$WORK_INCOME_TCR = as.factor(kidney_data$WORK_INCOME_TCR)
kidney_data$MED_COND_TRR = as.factor(kidney_data$MED_COND_TRR)
kidney_data$LOS = as.numeric(kidney_data$LOS)
kidney_data$TX_PROCEDUR_TY_KI = as.factor(kidney_data$TX_PROCEDUR_TY_KI)
kidney_data$INIT_STAT = as.factor(kidney_data$INIT_STAT)
kidney_data$END_STAT = as.factor(kidney_data$END_STAT)
kidney_data$HGT_CM_CALC = as.numeric(kidney_data$HGT_CM_CALC)
kidney_data$B1 <- as.factor(kidney_data$B1)
kidney_data$DB1 <- as.factor(kidney_data$DB1)
kidney_data$B2 <- as.factor(kidney_data$B2)
kidney_data$DB2 <- as.factor(kidney_data$DB2)
kidney_data$DDR1 <- as.factor(kidney_data$DDR1)
kidney_data$DR1 <- as.factor(kidney_data$DR1)
kidney_data$DDR1 <- as.factor(kidney_data$DDR1)
kidney_data$DR2 <- as.factor(kidney_data$DR2)
kidney_data$DDR2 <- as.factor(kidney_data$DDR2)
kidney_data$ABO_incompatible <- as.factor(kidney_data$ABO_incompatible)

kidney_data <- within(kidney_data, {
  ABO_A <- as.numeric(ABO == "A")
  ABO_B <- as.numeric(ABO == "B")
  ABO_AB <- as.numeric(ABO == "AB")
  ABO_O <- as.numeric(ABO == "O") 
})


kidney_data <- kidney_data %>%
  mutate(GRF_FAIL_CAUSE_TY_KI = case_when(
    GRF_FAIL_CAUSE_TY_KI == 1 ~ "Hyperacute Rejection",
    GRF_FAIL_CAUSE_TY_KI == 2 ~ "Acute Rejection",
    GRF_FAIL_CAUSE_TY_KI == 3 ~ "Primary Failure",
    GRF_FAIL_CAUSE_TY_KI == 4 ~ "Graft Thrombosis",
    GRF_FAIL_CAUSE_TY_KI == 5 ~ "Infection",
    GRF_FAIL_CAUSE_TY_KI == 6 ~ "Surgical Complications",
    GRF_FAIL_CAUSE_TY_KI == 7 ~ "Urological Complications",
    GRF_FAIL_CAUSE_TY_KI == 8 ~ "Recurrent Disease",
    GRF_FAIL_CAUSE_TY_KI == 9 ~ "Primary Non-Function (Graft Never Functioned Post-Transplant)", # Note: This is repeated for code 12
    GRF_FAIL_CAUSE_TY_KI == 10 ~ "Chronic Rejection",
    GRF_FAIL_CAUSE_TY_KI == 11 ~ "BK (Polyoma) Virus",
    GRF_FAIL_CAUSE_TY_KI == 12 ~ "Primary Non-Function (Graft Never Functioned Post-Transplant)",
    GRF_FAIL_CAUSE_TY_KI == 999 ~ "Other",
    TRUE ~ "Unknown" 
  ))

kidney_data$GRF_FAIL_CAUSE_TY_KI = as.factor(kidney_data$GRF_FAIL_CAUSE_TY_KI)

kidney_data <- kidney_data %>% 
  mutate(DWFG_KI = case_when(
    DWFG_KI == "Y" ~ 1,
    DWFG_KI == "N" ~ 0
  ))

kidney_data$DWFG_KI <- as.factor(kidney_data$DWFG_KI)
kidney_data$GSTATUS_KI <- as.numeric(kidney_data$GSTATUS_KI)
kidney_data <- kidney_data %>% filter(ABO %in% c("A", "B", "AB", "O"))
kidney_data$ABO <- droplevels(as.factor(kidney_data$ABO))
kidney_data$BMI_DON_CALC <- as.numeric(kidney_data$BMI_DON_CALC)
kidney_data$PERM_STATE <- as.factor(kidney_data$PERM_STATE)
kidney_data$CITIZENSHIP <- as.factor(kidney_data$CITIZENSHIP)
kidney_data$PRI_PAYMENT_TCR_KI <- as.factor(kidney_data$PRI_PAYMENT_TCR_KI)
kidney_data <- kidney_data %>%
  mutate(PRI_PAYMENT_TCR_KI = case_when(
    PRI_PAYMENT_TCR_KI %in% c(3, 4, 13) ~ "Medicare",
    PRI_PAYMENT_TCR_KI == 1 ~ "Private Insurance",
    PRI_PAYMENT_TCR_KI == 2 ~ "Medicaid",
    PRI_PAYMENT_TCR_KI == 5 ~ "CHIP",
    PRI_PAYMENT_TCR_KI == 6 ~ "VA",
    PRI_PAYMENT_TCR_KI == 7 ~ "Other Gov",
    PRI_PAYMENT_TCR_KI == 8 ~ "Self",
    PRI_PAYMENT_TCR_KI %in% c(9, 10, 11, 12) ~ "Other",
    TRUE ~ NA_character_
  ))
kidney_data$PRI_PAYMENT_TCR_KI <- as.factor(kidney_data$PRI_PAYMENT_TCR_KI)

kidney_data <- kidney_data %>%
  mutate(CITIZENSHIP = case_when(
    CITIZENSHIP == 1 ~ "US Citizen",
    CITIZENSHIP == 2 ~ "RESIDENT ALIEN",
    CITIZENSHIP %in% c(3,5,6) ~ "NON-RESIDENT ALIEN",
    CITIZENSHIP == 4 ~ "Non-Citizen, US Resident",
    TRUE ~ NA_character_  
  ))

kidney_data$CITIZENSHIP <- as.factor(kidney_data$CITIZENSHIP)
kidney_data$LISTING_CTR_CODE <- as.factor(kidney_data$LISTING_CTR_CODE)
kidney_data$COD_KI <- as.factor(kidney_data$COD_KI)
kidney_data <- kidney_data %>% mutate(COD_KI = case_when(
  COD_KI %in% c("3200", "3201", "3202", "3203" , "3204", "3299") ~ "Graft Fail", 
  COD_KI %in% c("3300", "3301", "3302", "3303", "3304", "3305", "3306", "3307", "3308", "3399") ~ "Infection",
  COD_KI %in% c("3400", "3401", "3402", "3499") ~ "Cardiovascular",
  COD_KI %in% c("3500", "3599") ~ "Cerebrovascular",
  COD_KI %in% c("3600", "3601", "3699") ~ "Hemorrhage",
  COD_KI %in% c("3700","3701","3702","3799") ~ "Malignancy",
  COD_KI %in% c("3800","3899") ~ "Trauma",
  COD_KI %in% c("3900","3901","3902","3903", "3904", "3905", "3906", "3907","3908", "3909", "3910", "3911", "3912",
                "3913","3914") ~ "Misc.",
  COD_KI == "3915" ~ "Primary non-function",
  COD_KI %in% c("3916","3917") ~ "Viral Infection",
  COD_KI == "."  ~ "Unknown"
))

kidney_data$CTR_CODE <- as.factor(kidney_data$CTR_CODE)

kidney_data <- kidney_data %>%
  mutate(
    TRTREJ1Y_KI = fct_na_value_to_level(as.factor(TRTREJ1Y_KI), level = "U"),
    TRTREJ6M_KI = fct_na_value_to_level(as.factor(TRTREJ6M_KI), level = "U"),
    WORK_INCOME_TCR = fct_na_value_to_level(as.factor(WORK_INCOME_TCR), level = "U"),
    WORK_INCOME_TRR = fct_na_value_to_level(as.factor(WORK_INCOME_TRR), level = "U")
  )



kidney_data$FUNC_STAT_TRF <- as.factor(kidney_data$FUNC_STAT_TRF)


dmu_data <- data.frame(
  Group = kidney_data$ETHCAT,
  ETHCAT = kidney_data$ETHCAT,
  WaitlistDuration = kidney_data$DAYSWAIT_CHRON_KI, # Input 1
  QualityScore = kidney_data$LKDPI, # Input 2
  OutcomeScore = kidney_data$GTIME_KI, # Output 1
  WL_days = kidney_data$DAYSWAIT_CHRON_KI,
  LKDPI = kidney_data$LKDPI,
  GTIME_KI = kidney_data$GTIME_KI,
  Year = kidney_data$TX_YEAR,
  DIABETES_DON = kidney_data$DIABETES_DON,
  eGFR = kidney_data$eGFR,
  BMI_CALC = kidney_data$BMI,
  WGT_KG_CALC = kidney_data$WGT_KG_CALC, 
  ABO = kidney_data$ABO, 
  ON_DIALYSIS = kidney_data$ON_DIALYSIS,
  AGE = kidney_data$AGE, 
  AGE_DON = kidney_data$AGE_DON, 
  HGT_CM_DON_CALC =  kidney_data$HGT_CM_DON_CALC, 
  WGT_KG_DON_CALC = kidney_data$WGT_KG_DON_CALC, 
  ETHCAT_DON = kidney_data$ETHCAT_DON, 
  HIST_HYPERTENS_DON = kidney_data$HIST_HYPERTENS_DON,
  PRA = kidney_data$PRA, 
  GENDER = kidney_data$GENDER, 
  DISTANCE = kidney_data$DISTANCE, 
  PRI_PAYMENT_TRR_KI = kidney_data$PRI_PAYMENT_TRR_KI, 
  PREV_KI_TX = kidney_data$PREV_KI_TX, 
  HIST_DIABETES_DON = kidney_data$HIST_DIABETES_DON,
  CTR_CODE = kidney_data$CTR_CODE,
  DIAG_KI = kidney_data$DIAG_KI, 
  BMI_DON_CALC = kidney_data$BMI_DON_CALC,
  GRF_FAIL_CAUSE_TY_KI = kidney_data$GRF_FAIL_CAUSE_TY_KI,
  ACADEMIC_LEVEL_TCR = kidney_data$ACADEMIC_LEVEL_TCR,
  DIAB = kidney_data$DIAB,
  EDUCATION = kidney_data$EDUCATION,
  ABO_MAT = kidney_data$ABO_MAT,
  PREV_TX_ANY = kidney_data$PREV_TX_ANY, 
  REGION = kidney_data$REGION, 
  MED_COND_TRR = kidney_data$MED_COND_TRR, 
  LOS = kidney_data$LOS,
  INIT_STAT = kidney_data$INIT_STAT, 
  END_STAT = kidney_data$END_STAT, 
  NPKID = kidney_data$NPKID,
  ABO_A = kidney_data$ABO_A,
  ABO_B = kidney_data$ABO_B,
  ABO_AB = kidney_data$ABO_B,
  ABO_O = kidney_data$ABO_O, 
  ETH_Asian = kidney_data$ETH_Asian,
  ETH_Black = kidney_data$ETH_Black,
  ETH_Hispanic = kidney_data$ETH_Hispanic,
  ETH_White = kidney_data$ETH_White,
  DWFG_KI = kidney_data$DWFG_KI,
  GSTATUS_KI = kidney_data$GSTATUS_KI,
  HIST_CIG_DON = kidney_data$HIST_CIG_DON,
  SBP = kidney_data$SBP,
  GENDER_DON = kidney_data$GENDER_DON,
  B1 = kidney_data$B1,
  DB1 = kidney_data$DB1,
  B2 = kidney_data$B2,
  DB2 = kidney_data$DB2,
  DDR1 = kidney_data$DDR1,
  DR1 = kidney_data$DR1,
  DDR2 = kidney_data$DDR2,
  DR2 = kidney_data$DR2,
  HLA_DR_mismatches = kidney_data$HLA_DR_mismatches,
  HLA_B_mismatches = kidney_data$HLA_B_mismatches,
  ABO_incompatible = kidney_data$ABO_incompatible,
  EDUCATION = kidney_data$EDUCATION,
  CITIZENSHIP = kidney_data$CITIZENSHIP,
  PRI_PAYMENT_TCR_KI = kidney_data$PRI_PAYMENT_TCR_KI,
  LISTING_CTR_CODE = kidney_data$LISTING_CTR_CODE,
  WORK_INCOME_TRR = kidney_data$WORK_INCOME_TRR,
  WORK_INCOME_TCR = kidney_data$WORK_INCOME_TCR,
  TRTREJ1Y_KI = kidney_data$TRTREJ1Y_KI, 
  TRTREJ6M_KI = kidney_data$TRTREJ6M_KI,
  FUNC_STAT_TRF = kidney_data$FUNC_STAT_TRF,
  PREV_KI_TX = kidney_data$PREV_KI_TX
)

dmu_data <- dmu_data[complete.cases(dmu_data),]

dim(dmu_data)
table(dmu_data$Group)/nrow(dmu_data)

write.csv(dmu_data, "dmu_data.csv", row.names = F) # For Simulation in Supplementary Materials (S2) 

####################################################################################################
#########################################  RESAMPLING (Section 3.2) ################################
####################################################################################################

# ESRD Prevalence Counts by Year (https://usrds-adr.niddk.nih.gov/2024/reference-tables, Table B.1)
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

# Convert counts to proportions
get_proportions <- function(counts) {
  total <- sum(counts)
  round(counts / total, 5)
}

# Stratified sampling function
stratified_sample <- function(data, group_var, desired_props, total_samples) {
  desired_counts <- round(desired_props[names(table(data[[group_var]]))] * total_samples)
  sampled_data <- map2_dfr(split(data, data[[group_var]]), names(desired_counts), function(group_data, group_name) {
    slice_sample(group_data, n = min(nrow(group_data), desired_counts[group_name]))
  })
  return(sampled_data)
}

set.seed(14)
# Perform resampling for each year
years <- intersect(names(esrd_counts_by_year), unique(dmu_data$Year) |> as.character())
balanced_data_list <- lapply(years, function(year) {
  year_data <- dmu_data %>% filter(Year == as.integer(year))
  if (nrow(year_data) == 0) return(NULL)
  
  desired_props <- get_proportions(esrd_counts_by_year[[year]])
  group_sizes <- table(year_data$Group)
  max_possible_samples <- min(group_sizes / desired_props[names(group_sizes)])
  total_samples <- floor(max_possible_samples)
  
  stratified_sample(year_data, "Group", desired_props, total_samples)
})

# Combine all years into one dataset
Final_Data <- bind_rows(balanced_data_list)


drop_all_levels <- function(data) {
  col_classes <- sapply(data, class)
  categorical_cols <- names(data)[col_classes %in% c("factor", "character")]
  for (col in categorical_cols) {
    if (is.character(data[[col]])) {
      data[[col]] <- factor(data[[col]])
    }
    data[[col]] <- droplevels(data[[col]])
  }
  return(data)
}

Final_Data <- drop_all_levels(Final_Data)


################################################################################
Final_Data <- Final_Data %>% group_by(Year) %>% 
  mutate(Priority_Score = WaitlistDuration - mean(WaitlistDuration),
         Access_Score = LKDPI - mean(LKDPI),
         Outcome_Score = GTIME_KI - mean(GTIME_KI)
  ) %>% ungroup()

write.csv(Final_Data, "Final_Data.csv", row.names = F)


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
Final_Data %>%
  group_by(Group) %>%
  summarise(
    Mean = formatC(mean(WaitlistDuration, na.rm = TRUE), format = "f", digits = 1),
    SD = formatC(sd(WaitlistDuration, na.rm = TRUE), format = "f", digits = 1),
    Median = formatC(median(WaitlistDuration, na.rm = TRUE), format = "f", digits = 2),
    Q1 = formatC(quantile(WaitlistDuration, 0.25, na.rm = TRUE), format = "f", digits = 2),
    Q3 = formatC(quantile(WaitlistDuration, 0.75, na.rm = TRUE), format = "f", digits = 2)
  )

# Raw LKDPI Summaries
Final_Data %>%
  group_by(Group) %>%
  summarise(
    Mean = formatC(mean(LKDPI, na.rm = TRUE), format = "f", digits = 2),
    SD = formatC(sd(LKDPI, na.rm = TRUE), format = "f", digits = 1),
    Median = formatC(median(LKDPI, na.rm = TRUE), format = "f", digits = 2),
    Q1 = formatC(quantile(LKDPI, 0.25, na.rm = TRUE), format = "f", digits = 2),
    Q3 = formatC(quantile(LKDPI, 0.75, na.rm = TRUE), format = "f", digits = 2)
  )


# Raw Graft Lifespan Summaries
Final_Data %>%
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
Final_Data %>%
  group_by(Group) %>%
  summarise(
    Mean = formatC(mean(Priority_Score, na.rm = TRUE), format = "f", digits = 1),
    SD = formatC(sd(Priority_Score, na.rm = TRUE), format = "f", digits = 1),
    Median = formatC(median(Priority_Score, na.rm = TRUE), format = "f", digits = 2),
    Q1 = formatC(quantile(Priority_Score, 0.25, na.rm = TRUE), format = "f", digits = 2),
    Q3 = formatC(quantile(Priority_Score, 0.75, na.rm = TRUE), format = "f", digits = 2)
  )

# Relative LKDPI Summaries
Final_Data %>%
  group_by(Group) %>%
  summarise(
    Mean = formatC(mean(Access_Score, na.rm = TRUE), format = "f", digits = 1),
    SD = formatC(sd(Access_Score, na.rm = TRUE), format = "f", digits = 1),
    Median = formatC(median(Access_Score, na.rm = TRUE), format = "f", digits = 2),
    Q1 = formatC(quantile(Access_Score, 0.25, na.rm = TRUE), format = "f", digits = 2),
    Q3 = formatC(quantile(Access_Score, 0.75, na.rm = TRUE), format = "f", digits = 2)
  )


# Relative Graft Lifespan Summaries
Final_Data %>%
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
Final_Data %>%
  group_by(Year, Group) %>%
  summarise(Average_Score = mean(Priority_Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(
    names_from = Year,
    values_from = Average_Score
  )


Final_Data %>%
  group_by(Year, Group) %>%
  summarise(Average_Score = mean(Access_Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(
    names_from = Year,
    values_from = Average_Score
  )

Final_Data %>%
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

# Sample Figure for 2019
Final_Data$Group <- factor(
  Final_Data$Group, levels = c("Asian", "Black", "Hispanic", "White")
)

# Prepare data
original_data <- dmu_data %>% 
  filter(Year == 2019) %>% 
  count(Group) %>%
  mutate(Percentage = n / sum(n),
         Dataset = "Imbalanced Data")

resampled_data <- Final_Data %>%
  filter(Year == 2019) %>% 
  count(Group) %>%
  mutate(Percentage = case_when(
    Group == "Asian" ~ 0.054,
    Group == "Black" ~ 0.300,
    Group == "Hispanic" ~ 0.191,
    Group == "White" ~ 0.455
  ),
  Dataset = "ESRD Prevalence")


combined_data <- rbind(original_data, resampled_data)
combined_data$Dataset <- factor(
  combined_data$Dataset,
  levels = c("Imbalanced Data", "ESRD Prevalence")
)


color_palette <- c("Asian" = "#08306b", "Black" = "#1f78b4",
                   "Hispanic" = "#6baed6", "White" = "#b3cde3")

ggplot(combined_data, aes(x = Dataset, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage * 100)), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 4, fontface = "bold") +  
  scale_fill_manual(values = color_palette) +
  labs(x = "", y = "Relative Frequency") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0, 0.75)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold"),  
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20, face = "bold", margin = ggplot2::margin(r = 10)),  
        legend.title = element_text(size = 18, face = "bold"), 
        legend.text = element_text(size = 16)) 
