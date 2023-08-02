# Get normalised and scaled gene expression matrix using cBioPortalData
library(cBioPortalData)
library(tidyverse)
library(reshape2)
library(survival)
library(survminer)

cohort <- "skcm_tcga" # TCGA cohort (skcm - skin cutaneous melanoma)
GOI <- c("AXL", "MITF") # the list genes of interests
GOI.survival = "AXL" # The gene to be tested for survival analysis, must be included in GOI

# Fetch expression data (normalised z-score) of your GOI
cbio <- cBioPortal()
skcm.exp.data <- getDataByGenes(cbio, 
                                studyId = cohort, 
                                molecularProfileIds = "skcm_tcga_rna_seq_v2_mrna_median_Zscores", 
                                genes = GOI, 
                                by = "hugoGeneSymbol")

# Change data to "wide" format
skcm.exp.data.wide <- skcm.exp.data$skcm_tcga_rna_seq_v2_mrna_median_Zscores %>%
  dcast(sampleId ~ hugoGeneSymbol, value.var = "value")

# Fetch clinical data
skcm_clinical <- clinicalData(cbio, "skcm_tcga") %>%
  left_join(skcm.exp.data.wide, by = "sampleId") %>% # join the clinical data and expression data
  mutate(Stage = ifelse(grepl("IV", AJCC_PATHOLOGIC_TUMOR_STAGE), "IV", # clean up stage data
                        ifelse(grepl("III", AJCC_PATHOLOGIC_TUMOR_STAGE), "III", 
                               ifelse(grepl("II", AJCC_PATHOLOGIC_TUMOR_STAGE), "II",
                                      ifelse(grepl("I", AJCC_PATHOLOGIC_TUMOR_STAGE), "I", NA))))) %>%
  mutate(Stage = factor(Stage, levels = c("I", "II", "III", "IV")))%>%
  dplyr::filter(OS_STATUS != "" & !is.na(OS_MONTHS)) %>% # filter for patients with OS data
  mutate(OS_STATUS_level = ifelse(OS_STATUS=="1:DECEASED",1,0)) 

# stratify patients by high & low GOI expression according to median
mrna_level = paste(GOI.survival, "level", sep = ".")
dataset <- skcm_clinical %>%
  dplyr::filter(!!rlang::sym(GOI.survival) != "NaN") %>% # filter out missing mrna expression level cases
  mutate(!!mrna_level := ifelse(!!rlang::sym(GOI.survival) > median(!!rlang::sym(GOI.survival)), 
                                "High (> median)",
                                "Low (< median)")) # Divide rna-expression level into two groups

# perform KM survival analysis
survival <- Surv(time = dataset$OS_MONTHS, event = dataset$OS_STATUS_level) 
frm <- paste("survival ~ ", mrna_level, sep = "") # Create a survival object and fit formula
km_mrna_fit <- survfit(as.formula(frm), data=dataset) # Fit the Kaplan-Meier curves

# Plot KM survival curve
ggsurv <- ggsurvplot(km_mrna_fit, data = dataset, pval = TRUE, pval.method = TRUE, font.x = c(18), font.y = c(18)) +
  labs(x = "Overall survival time (months)", colour = "")
ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 13, color = "black", face = "bold"), axis.title =  element_text(size = 20))
print(ggsurv$plot)

# scatter plot + calculate correlations of two genes
ggscatter(skcm_clinical, x = "MITF", y = "AXL", 
          conf.int = TRUE, 
          add.params = list(color="blue", fill="lightgrey")) +
  stat_cor(method = "spearman", label.x = 0.5) +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 8))

