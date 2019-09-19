
#### Read in proteinGroups and pre-process data ####
proteinGroups <- read.delim("proteinGroups.txt", stringsAsFactors=FALSE)

# Remove contaminants, reverse and only identified by site
proteinGroups <- proteinGroups[proteinGroups$Only.identified.by.site == "" ,]
proteinGroups <- proteinGroups[proteinGroups$Potential.contaminant == "" ,]
proteinGroups <- proteinGroups[proteinGroups$Reverse == "" ,]

# Process ratios
Ratios <- proteinGroups[, grep("Ratio.H.L.normalized.", colnames(proteinGroups))]
Ratios <- log2(Ratios)
median_Ratios <- apply(Ratios, 2, median, na.rm=TRUE)
Ratios <- sweep(Ratios, 2, median_Ratios, FUN="-")

# Other data subsets
Annotation <- proteinGroups[, c("Majority.protein.IDs", "Protein.names", "Gene.names")]
Intensities <- proteinGroups[, grep("Intensity...+", colnames(proteinGroups))]        # Get all columns reporting intensities
  Intensities <- Intensities[, -grep("Intensity\\.L\\..+", colnames(Intensities))]    # Drop reports of "light" intensities
  Intensities <- Intensities[, -grep("Intensity\\.H\\..+", colnames(Intensities))]    # Drop reports of "heavy" intensities
Ratio_counts <- proteinGroups[, grep("Ratio.H.L.count.", colnames(proteinGroups))]


#### Calculate significant, reproducible outliers ####

# Write out txt file for Perseus (and use it in Perseus to calculate significance B)
Perseus_data <- data.frame(proteinGroups$Majority.protein.IDs, Ratios, Intensities)
write.table(Perseus_data, "PerseusIn.txt", row.names = FALSE, sep="\t")

# Read in significanceB
PerseusOut <- read.delim("PerseusOut.txt", comment.char="#", stringsAsFactors=FALSE)
SigB <- PerseusOut[, grep("Significance.B", colnames(PerseusOut))]

# Is it a significant (P < 0.05) and reproducible (in at least 2 exp) outlier?
ATMi_A_sig.rep.outlier <- apply( SigB[, grep("Ratio.H.L.normalized.ATMi_A", colnames(SigB))], 1, function(x){ sum( x<0.05 ) >= 2 })  
ATMi_B_sig.rep.outlier <- apply( SigB[, grep("Ratio.H.L.normalized.ATMi_B", colnames(SigB))], 1, function(x){ sum( x<0.05 ) >= 2 })  
CHX_sig.rep.outlier <- apply( SigB[, grep("Ratio.H.L.normalized.CHX", colnames(SigB))], 1, function(x){ sum( x<0.05 ) >= 2 })
CPT_sig.rep.outlier <- apply( SigB[, grep("Ratio.H.L.normalized.CPT_", colnames(SigB))], 1, function(x){ sum( x<0.05 ) >= 2 })
HU_sig.rep.outlier <- apply( SigB[, grep("Ratio.H.L.normalized.HU", colnames(SigB))], 1, function(x){ sum( x<0.05 ) >= 2 })
NCC_sig.rep.outlier <- apply( SigB[, grep("Ratio.H.L.normalized.NCC", colnames(SigB))], 1, function(x){ sum( x<0.05 ) >= 2 })
Rosco_sig.rep.outlier <- apply( SigB[, grep("Ratio.H.L.normalized.Rosco_NCC", colnames(SigB))], 1, function(x){ sum( x<0.05 ) >= 2 })
RoscoCPT_sig.rep.outlier <- apply( SigB[, grep("Ratio.H.L.normalized.RoscoCPT", colnames(SigB))], 1, function(x){ sum( x<0.05 ) >= 2 })
TSA_sig.rep.outlier <- apply( SigB[, grep("Ratio.H.L.normalized.TSA", colnames(SigB))], 1, function(x){ sum( x<0.05 ) >= 2 })

Significance <- data.frame(PerseusOut$proteinGroups.Majority.protein.IDs, ATMi_A_sig.rep.outlier, ATMi_B_sig.rep.outlier,
                           CHX_sig.rep.outlier, CPT_sig.rep.outlier, HU_sig.rep.outlier, NCC_sig.rep.outlier, Rosco_sig.rep.outlier, 
                           RoscoCPT_sig.rep.outlier, TSA_sig.rep.outlier)

#### Assemble result file ####
Results <- data.frame(Annotation, Ratios, Ratio_counts, Intensities)

# Append significance 
Results <- merge(Results, Significance, by.x="Majority.protein.IDs", by.y="PerseusOut.proteinGroups.Majority.protein.IDs")

rm( list = ls()[! ls() %in% c("Results")] )   # Clear workspace


#### Add ICPs ####
library(readxl)
ICPs <- read_excel("Fuzzy_table_S1.xlsx", sheet = 1, col_names = TRUE)
ICPs$SimpleID <- gsub("-.+", "", ICPs$`Majority protein IDs`)                  # Simplify IDs
ICPs$SimpleID <- gsub(";.+", "", ICPs$SimpleID) 
ICPs <- ICPs[ ICPs$SimpleID %in% names(which(table(ICPs$SimpleID) == 1)) ,]    # Remove duplicates
ICPs <- ICPs[, c(10,2,9) ]

# Correct a wrong spelling
ICPs$`Lab internal category` <- gsub("chromatin structure & organization", "Chromatin structure & organization", ICPs$`Lab internal category`)

Results$SimpleID <- gsub("-.+", "", Results$Majority.protein.IDs)              # Create similar simplified IDs for the Results df
Results$SimpleID <- gsub(";.+", "", Results$SimpleID)

Results <- merge(Results, ICPs, all.x=TRUE)


#### Clean-up column names ####

# Unify replicate naming
colnames(Results) <- gsub("CHX([[:digit:]])", "CHX_\\1", colnames(Results))
colnames(Results) <- gsub("CPT_s([[:digit:]])", "CPT_\\1", colnames(Results))
colnames(Results) <- gsub("HU_Ex([[:digit:]])", "HU_\\1", colnames(Results))
colnames(Results) <- gsub("NCC_A", "NCC_1" , colnames(Results))
colnames(Results) <- gsub("NCC_B", "NCC_2" , colnames(Results))
colnames(Results) <- gsub("NCC_C", "NCC_3" , colnames(Results))
colnames(Results) <- gsub("Rosco_NCC36", "Rosco_1" , colnames(Results))
colnames(Results) <- gsub("Rosco_NCC37", "Rosco_2" , colnames(Results))
colnames(Results) <- gsub("Rosco_NCC38", "Rosco_3" , colnames(Results))
colnames(Results) <- gsub("TSA_Exp([[:digit:]])", "TSA_\\1", colnames(Results))

# More clean-up
colnames(Results) <- gsub(" ", "_", colnames(Results))
colnames(Results) <- gsub("Ratio.H.L.normalized.", "log2.SILAC.ratio.", colnames(Results))
colnames(Results) <- gsub("Ratio.H.L.count.", "Ratio.count.", colnames(Results))
colnames(Results) <- gsub(".", "_", colnames(Results), fixed = TRUE)



#### Calculate mean and SEM over triplicates ####

## For ATMi_A experiment ##
# For ratios with >= 1 ratio count
R1 <- ifelse( Results$Ratio_count_ATMi_A_1 >= 1 , Results$log2_SILAC_ratio_ATMi_A_1 , NA  )
R2 <- ifelse( Results$Ratio_count_ATMi_A_2 >= 1 , Results$log2_SILAC_ratio_ATMi_A_2 , NA  )
R3 <- ifelse( Results$Ratio_count_ATMi_A_3 >= 1 , Results$log2_SILAC_ratio_ATMi_A_3 , NA  )
ATMi_A_mean_ratio_1rc <- apply( data.frame(R1, R2, R3) , 1, mean)
ATMi_A_SEM_1rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

# For ratios with >= 2 ratio count
R1 <- ifelse( Results$Ratio_count_ATMi_A_1 >= 2 , Results$log2_SILAC_ratio_ATMi_A_1 , NA  )
R2 <- ifelse( Results$Ratio_count_ATMi_A_2 >= 2 , Results$log2_SILAC_ratio_ATMi_A_2 , NA  )
R3 <- ifelse( Results$Ratio_count_ATMi_A_3 >= 2 , Results$log2_SILAC_ratio_ATMi_A_3 , NA  )
ATMi_A_mean_ratio_2rc <- apply( data.frame(R1, R2, R3) , 1, mean)
ATMi_A_SEM_2rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

R1 <- ifelse( Results$Ratio_count_ATMi_A_1 >= 3 , Results$log2_SILAC_ratio_ATMi_A_1 , NA  )
R2 <- ifelse( Results$Ratio_count_ATMi_A_2 >= 3 , Results$log2_SILAC_ratio_ATMi_A_2 , NA  )
R3 <- ifelse( Results$Ratio_count_ATMi_A_3 >= 3 , Results$log2_SILAC_ratio_ATMi_A_3 , NA  )
ATMi_A_mean_ratio_3rc <- apply( data.frame(R1, R2, R3) , 1, mean)
ATMi_A_SEM_3rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)


## For ATMi_B experiment ##
# For ratios with >= 1 ratio count
R1 <- ifelse( Results$Ratio_count_ATMi_B_1 >= 1 , Results$log2_SILAC_ratio_ATMi_B_1 , NA  )
R2 <- ifelse( Results$Ratio_count_ATMi_B_2 >= 1 , Results$log2_SILAC_ratio_ATMi_B_2 , NA  )
R3 <- ifelse( Results$Ratio_count_ATMi_B_3 >= 1 , Results$log2_SILAC_ratio_ATMi_B_3 , NA  )
ATMi_B_mean_ratio_1rc <- apply( data.frame(R1, R2, R3) , 1, mean)
ATMi_B_SEM_1rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

# For ratios with >= 2 ratio count
R1 <- ifelse( Results$Ratio_count_ATMi_B_1 >= 2 , Results$log2_SILAC_ratio_ATMi_B_1 , NA  )
R2 <- ifelse( Results$Ratio_count_ATMi_B_2 >= 2 , Results$log2_SILAC_ratio_ATMi_B_2 , NA  )
R3 <- ifelse( Results$Ratio_count_ATMi_B_3 >= 2 , Results$log2_SILAC_ratio_ATMi_B_3 , NA  )
ATMi_B_mean_ratio_2rc <- apply( data.frame(R1, R2, R3) , 1, mean)
ATMi_B_SEM_2rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

R1 <- ifelse( Results$Ratio_count_ATMi_B_1 >= 3 , Results$log2_SILAC_ratio_ATMi_B_1 , NA  )
R2 <- ifelse( Results$Ratio_count_ATMi_B_2 >= 3 , Results$log2_SILAC_ratio_ATMi_B_2 , NA  )
R3 <- ifelse( Results$Ratio_count_ATMi_B_3 >= 3 , Results$log2_SILAC_ratio_ATMi_B_3 , NA  )
ATMi_B_mean_ratio_3rc <- apply( data.frame(R1, R2, R3) , 1, mean)
ATMi_B_SEM_3rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)


## For CHX experiment ##
# For ratios with >= 1 ratio count
R1 <- ifelse( Results$Ratio_count_CHX_1 >= 1 , Results$log2_SILAC_ratio_CHX_1 , NA  )
R2 <- ifelse( Results$Ratio_count_CHX_2 >= 1 , Results$log2_SILAC_ratio_CHX_2 , NA  )
R3 <- ifelse( Results$Ratio_count_CHX_3 >= 1 , Results$log2_SILAC_ratio_CHX_3 , NA  )
CHX_mean_ratio_1rc <- apply( data.frame(R1, R2, R3) , 1, mean)
CHX_SEM_1rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

# For ratios with >= 2 ratio count
R1 <- ifelse( Results$Ratio_count_CHX_1 >= 2 , Results$log2_SILAC_ratio_CHX_1 , NA  )
R2 <- ifelse( Results$Ratio_count_CHX_2 >= 2 , Results$log2_SILAC_ratio_CHX_2 , NA  )
R3 <- ifelse( Results$Ratio_count_CHX_3 >= 2 , Results$log2_SILAC_ratio_CHX_3 , NA  )
CHX_mean_ratio_2rc <- apply( data.frame(R1, R2, R3) , 1, mean)
CHX_SEM_2rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

R1 <- ifelse( Results$Ratio_count_CHX_1 >= 3 , Results$log2_SILAC_ratio_CHX_1 , NA  )
R2 <- ifelse( Results$Ratio_count_CHX_2 >= 3 , Results$log2_SILAC_ratio_CHX_2 , NA  )
R3 <- ifelse( Results$Ratio_count_CHX_3 >= 3 , Results$log2_SILAC_ratio_CHX_3 , NA  )
CHX_mean_ratio_3rc <- apply( data.frame(R1, R2, R3) , 1, mean)
CHX_SEM_3rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)


## For CPT experiment ##
# For ratios with >= 1 ratio count
R1 <- ifelse( Results$Ratio_count_CPT_1 >= 1 , Results$log2_SILAC_ratio_CPT_1 , NA  )
R2 <- ifelse( Results$Ratio_count_CPT_2 >= 1 , Results$log2_SILAC_ratio_CPT_2 , NA  )
R3 <- ifelse( Results$Ratio_count_CPT_3 >= 1 , Results$log2_SILAC_ratio_CPT_3 , NA  )
CPT_mean_ratio_1rc <- apply( data.frame(R1, R2, R3) , 1, mean)
CPT_SEM_1rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

# For ratios with >= 2 ratio count
R1 <- ifelse( Results$Ratio_count_CPT_1 >= 2 , Results$log2_SILAC_ratio_CPT_1 , NA  )
R2 <- ifelse( Results$Ratio_count_CPT_2 >= 2 , Results$log2_SILAC_ratio_CPT_2 , NA  )
R3 <- ifelse( Results$Ratio_count_CPT_3 >= 2 , Results$log2_SILAC_ratio_CPT_3 , NA  )
CPT_mean_ratio_2rc <- apply( data.frame(R1, R2, R3) , 1, mean)
CPT_SEM_2rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

R1 <- ifelse( Results$Ratio_count_CPT_1 >= 3 , Results$log2_SILAC_ratio_CPT_1 , NA  )
R2 <- ifelse( Results$Ratio_count_CPT_2 >= 3 , Results$log2_SILAC_ratio_CPT_2 , NA  )
R3 <- ifelse( Results$Ratio_count_CPT_3 >= 3 , Results$log2_SILAC_ratio_CPT_3 , NA  )
CPT_mean_ratio_3rc <- apply( data.frame(R1, R2, R3) , 1, mean)
CPT_SEM_3rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)


## For HU experiment ##
# For ratios with >= 1 ratio count
R1 <- ifelse( Results$Ratio_count_HU_1 >= 1 , Results$log2_SILAC_ratio_HU_1 , NA  )
R2 <- ifelse( Results$Ratio_count_HU_2 >= 1 , Results$log2_SILAC_ratio_HU_2 , NA  )
R3 <- ifelse( Results$Ratio_count_HU_3 >= 1 , Results$log2_SILAC_ratio_HU_3 , NA  )
HU_mean_ratio_1rc <- apply( data.frame(R1, R2, R3) , 1, mean)
HU_SEM_1rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

# For ratios with >= 2 ratio count
R1 <- ifelse( Results$Ratio_count_HU_1 >= 2 , Results$log2_SILAC_ratio_HU_1 , NA  )
R2 <- ifelse( Results$Ratio_count_HU_2 >= 2 , Results$log2_SILAC_ratio_HU_2 , NA  )
R3 <- ifelse( Results$Ratio_count_HU_3 >= 2 , Results$log2_SILAC_ratio_HU_3 , NA  )
HU_mean_ratio_2rc <- apply( data.frame(R1, R2, R3) , 1, mean)
HU_SEM_2rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

R1 <- ifelse( Results$Ratio_count_HU_1 >= 3 , Results$log2_SILAC_ratio_HU_1 , NA  )
R2 <- ifelse( Results$Ratio_count_HU_2 >= 3 , Results$log2_SILAC_ratio_HU_2 , NA  )
R3 <- ifelse( Results$Ratio_count_HU_3 >= 3 , Results$log2_SILAC_ratio_HU_3 , NA  )
HU_mean_ratio_3rc <- apply( data.frame(R1, R2, R3) , 1, mean)
HU_SEM_3rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)


## For NCC experiment ##
# For ratios with >= 1 ratio count
R1 <- ifelse( Results$Ratio_count_NCC_1 >= 1 , Results$log2_SILAC_ratio_NCC_1 , NA  )
R2 <- ifelse( Results$Ratio_count_NCC_2 >= 1 , Results$log2_SILAC_ratio_NCC_2 , NA  )
R3 <- ifelse( Results$Ratio_count_NCC_3 >= 1 , Results$log2_SILAC_ratio_NCC_3 , NA  )
NCC_mean_ratio_1rc <- apply( data.frame(R1, R2, R3) , 1, mean)
NCC_SEM_1rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

# For ratios with >= 2 ratio count
R1 <- ifelse( Results$Ratio_count_NCC_1 >= 2 , Results$log2_SILAC_ratio_NCC_1 , NA  )
R2 <- ifelse( Results$Ratio_count_NCC_2 >= 2 , Results$log2_SILAC_ratio_NCC_2 , NA  )
R3 <- ifelse( Results$Ratio_count_NCC_3 >= 2 , Results$log2_SILAC_ratio_NCC_3 , NA  )
NCC_mean_ratio_2rc <- apply( data.frame(R1, R2, R3) , 1, mean)
NCC_SEM_2rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

R1 <- ifelse( Results$Ratio_count_NCC_1 >= 3 , Results$log2_SILAC_ratio_NCC_1 , NA  )
R2 <- ifelse( Results$Ratio_count_NCC_2 >= 3 , Results$log2_SILAC_ratio_NCC_2 , NA  )
R3 <- ifelse( Results$Ratio_count_NCC_3 >= 3 , Results$log2_SILAC_ratio_NCC_3 , NA  )
NCC_mean_ratio_3rc <- apply( data.frame(R1, R2, R3) , 1, mean)
NCC_SEM_3rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)


## For Rosco experiment ##
# For ratios with >= 1 ratio count
R1 <- ifelse( Results$Ratio_count_Rosco_1 >= 1 , Results$log2_SILAC_ratio_Rosco_1 , NA  )
R2 <- ifelse( Results$Ratio_count_Rosco_2 >= 1 , Results$log2_SILAC_ratio_Rosco_2 , NA  )
R3 <- ifelse( Results$Ratio_count_Rosco_3 >= 1 , Results$log2_SILAC_ratio_Rosco_3 , NA  )
Rosco_mean_ratio_1rc <- apply( data.frame(R1, R2, R3) , 1, mean)
Rosco_SEM_1rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

# For ratios with >= 2 ratio count
R1 <- ifelse( Results$Ratio_count_Rosco_1 >= 2 , Results$log2_SILAC_ratio_Rosco_1 , NA  )
R2 <- ifelse( Results$Ratio_count_Rosco_2 >= 2 , Results$log2_SILAC_ratio_Rosco_2 , NA  )
R3 <- ifelse( Results$Ratio_count_Rosco_3 >= 2 , Results$log2_SILAC_ratio_Rosco_3 , NA  )
Rosco_mean_ratio_2rc <- apply( data.frame(R1, R2, R3) , 1, mean)
Rosco_SEM_2rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

R1 <- ifelse( Results$Ratio_count_Rosco_1 >= 3 , Results$log2_SILAC_ratio_Rosco_1 , NA  )
R2 <- ifelse( Results$Ratio_count_Rosco_2 >= 3 , Results$log2_SILAC_ratio_Rosco_2 , NA  )
R3 <- ifelse( Results$Ratio_count_Rosco_3 >= 3 , Results$log2_SILAC_ratio_Rosco_3 , NA  )
Rosco_mean_ratio_3rc <- apply( data.frame(R1, R2, R3) , 1, mean)
Rosco_SEM_3rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)


## For RoscoCPT experiment ##
# For ratios with >= 1 ratio count
R1 <- ifelse( Results$Ratio_count_RoscoCPT_1 >= 1 , Results$log2_SILAC_ratio_RoscoCPT_1 , NA  )
R2 <- ifelse( Results$Ratio_count_RoscoCPT_2 >= 1 , Results$log2_SILAC_ratio_RoscoCPT_2 , NA  )
R3 <- ifelse( Results$Ratio_count_RoscoCPT_3 >= 1 , Results$log2_SILAC_ratio_RoscoCPT_3 , NA  )
RoscoCPT_mean_ratio_1rc <- apply( data.frame(R1, R2, R3) , 1, mean)
RoscoCPT_SEM_1rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

# For ratios with >= 2 ratio count
R1 <- ifelse( Results$Ratio_count_RoscoCPT_1 >= 2 , Results$log2_SILAC_ratio_RoscoCPT_1 , NA  )
R2 <- ifelse( Results$Ratio_count_RoscoCPT_2 >= 2 , Results$log2_SILAC_ratio_RoscoCPT_2 , NA  )
R3 <- ifelse( Results$Ratio_count_RoscoCPT_3 >= 2 , Results$log2_SILAC_ratio_RoscoCPT_3 , NA  )
RoscoCPT_mean_ratio_2rc <- apply( data.frame(R1, R2, R3) , 1, mean)
RoscoCPT_SEM_2rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

R1 <- ifelse( Results$Ratio_count_RoscoCPT_1 >= 3 , Results$log2_SILAC_ratio_RoscoCPT_1 , NA  )
R2 <- ifelse( Results$Ratio_count_RoscoCPT_2 >= 3 , Results$log2_SILAC_ratio_RoscoCPT_2 , NA  )
R3 <- ifelse( Results$Ratio_count_RoscoCPT_3 >= 3 , Results$log2_SILAC_ratio_RoscoCPT_3 , NA  )
RoscoCPT_mean_ratio_3rc <- apply( data.frame(R1, R2, R3) , 1, mean)
RoscoCPT_SEM_3rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)


## For TSA experiment ##
# For ratios with >= 1 ratio count
R1 <- ifelse( Results$Ratio_count_TSA_1 >= 1 , Results$log2_SILAC_ratio_TSA_1 , NA  )
R2 <- ifelse( Results$Ratio_count_TSA_2 >= 1 , Results$log2_SILAC_ratio_TSA_2 , NA  )
R3 <- ifelse( Results$Ratio_count_TSA_3 >= 1 , Results$log2_SILAC_ratio_TSA_3 , NA  )
TSA_mean_ratio_1rc <- apply( data.frame(R1, R2, R3) , 1, mean)
TSA_SEM_1rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

# For ratios with >= 2 ratio count
R1 <- ifelse( Results$Ratio_count_TSA_1 >= 2 , Results$log2_SILAC_ratio_TSA_1 , NA  )
R2 <- ifelse( Results$Ratio_count_TSA_2 >= 2 , Results$log2_SILAC_ratio_TSA_2 , NA  )
R3 <- ifelse( Results$Ratio_count_TSA_3 >= 2 , Results$log2_SILAC_ratio_TSA_3 , NA  )
TSA_mean_ratio_2rc <- apply( data.frame(R1, R2, R3) , 1, mean)
TSA_SEM_2rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)

R1 <- ifelse( Results$Ratio_count_TSA_1 >= 3 , Results$log2_SILAC_ratio_TSA_1 , NA  )
R2 <- ifelse( Results$Ratio_count_TSA_2 >= 3 , Results$log2_SILAC_ratio_TSA_2 , NA  )
R3 <- ifelse( Results$Ratio_count_TSA_3 >= 3 , Results$log2_SILAC_ratio_TSA_3 , NA  )
TSA_mean_ratio_3rc <- apply( data.frame(R1, R2, R3) , 1, mean)
TSA_SEM_3rc <-        apply( data.frame(R1, R2, R3) , 1, sd) / sqrt(3)


# Combine ratios and SEMs into a data frame
ratio_df <- data.frame( ATMi_A_mean_ratio_1rc, ATMi_A_SEM_1rc,
                        ATMi_A_mean_ratio_2rc, ATMi_A_SEM_2rc,
                        ATMi_A_mean_ratio_3rc, ATMi_A_SEM_3rc,
                      
                        ATMi_B_mean_ratio_1rc, ATMi_B_SEM_1rc,
                        ATMi_B_mean_ratio_2rc, ATMi_B_SEM_2rc,
                        ATMi_B_mean_ratio_3rc, ATMi_B_SEM_3rc,
                        
                        CHX_mean_ratio_1rc, CHX_SEM_1rc,
                        CHX_mean_ratio_2rc, CHX_SEM_2rc,
                        CHX_mean_ratio_3rc, CHX_SEM_3rc,
                        
                        CPT_mean_ratio_1rc, CPT_SEM_1rc,
                        CPT_mean_ratio_2rc, CPT_SEM_2rc,
                        CPT_mean_ratio_3rc, CPT_SEM_3rc,
                        
                        HU_mean_ratio_1rc, HU_SEM_1rc,
                        HU_mean_ratio_2rc, HU_SEM_2rc,
                        HU_mean_ratio_3rc, HU_SEM_3rc,
                        
                        NCC_mean_ratio_1rc, NCC_SEM_1rc,
                        NCC_mean_ratio_2rc, NCC_SEM_2rc,
                        NCC_mean_ratio_3rc, NCC_SEM_3rc,
                        
                        Rosco_mean_ratio_1rc, Rosco_SEM_1rc,
                        Rosco_mean_ratio_2rc, Rosco_SEM_2rc,
                        Rosco_mean_ratio_3rc, Rosco_SEM_3rc,
                        
                        RoscoCPT_mean_ratio_1rc, RoscoCPT_SEM_1rc,
                        RoscoCPT_mean_ratio_2rc, RoscoCPT_SEM_2rc,
                        RoscoCPT_mean_ratio_3rc, RoscoCPT_SEM_3rc,
                        
                        TSA_mean_ratio_1rc, TSA_SEM_1rc,
                        TSA_mean_ratio_2rc, TSA_SEM_2rc,
                        TSA_mean_ratio_3rc, TSA_SEM_3rc  )


#### Assemble final result file ####
# I removed the individual sums and added the individual intensities
Results <- data.frame(Results, ratio_df)


#### Write out results ####

# Complete results file
write.csv(Results, "NCC9b_results_Aug17.csv", row.names = FALSE)

