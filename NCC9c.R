library(data.table)

# Read in the last file
DT <- fread("NCC9b_results_Aug17.csv")

#### Add an "Replicate" column as in the first NCC paper ####

# This column should show in which replicates a protein had been detected with >= 2 ratio counts

# For ATMi_A
ATMi_A <- data.table( ifelse( DT$Ratio_count_ATMi_A_1 >= 2 , "1", "" ),   # If the ratio count is >= 2, return the replicate number
                      ifelse( DT$Ratio_count_ATMi_A_2 >= 2 , "2", "" ),
                      ifelse( DT$Ratio_count_ATMi_A_3 >= 2 , "3", "" ))
DT[, Replicate_2rc_ATMi_A := ATMi_A[, paste(V1, V2, V3, sep = "") ]]     # Paste replicate numbers into one string and append to the big table

# For ATMi_B
ATMi_B <- data.table( ifelse( DT$Ratio_count_ATMi_B_1 >= 2 , "1", "" ),   # If the ratio count is >= 2, return the replicate number
                      ifelse( DT$Ratio_count_ATMi_B_2 >= 2 , "2", "" ),
                      ifelse( DT$Ratio_count_ATMi_B_3 >= 2 , "3", "" ))
DT[, Replicate_2rc_ATMi_B := ATMi_B[, paste(V1, V2, V3, sep = "") ]]     # Paste replicate numbers into one string and append to the big table

# For CHX
CHX <- data.table( ifelse( DT$Ratio_count_CHX_1 >= 2 , "1", "" ),   # If the ratio count is >= 2, return the replicate number
                   ifelse( DT$Ratio_count_CHX_2 >= 2 , "2", "" ),
                   ifelse( DT$Ratio_count_CHX_3 >= 2 , "3", "" ))
DT[, Replicate_2rc_CHX := CHX[, paste(V1, V2, V3, sep = "") ]]     # Paste replicate numbers into one string and append to the big table

# For CPT
CPT <- data.table( ifelse( DT$Ratio_count_CPT_1 >= 2 , "1", "" ),   # If the ratio count is >= 2, return the replicate number
                   ifelse( DT$Ratio_count_CPT_2 >= 2 , "2", "" ),
                   ifelse( DT$Ratio_count_CPT_3 >= 2 , "3", "" ))
DT[, Replicate_2rc_CPT := CPT[, paste(V1, V2, V3, sep = "") ]]     # Paste replicate numbers into one string and append to the big table

# For HU
HU <- data.table( ifelse( DT$Ratio_count_HU_1 >= 2 , "1", "" ),   # If the ratio count is >= 2, return the replicate number
                  ifelse( DT$Ratio_count_HU_2 >= 2 , "2", "" ),
                  ifelse( DT$Ratio_count_HU_3 >= 2 , "3", "" ))
DT[, Replicate_2rc_HU := HU[, paste(V1, V2, V3, sep = "") ]]     # Paste replicate numbers into one string and append to the big table

# For NCC
NCC <- data.table( ifelse( DT$Ratio_count_NCC_1 >= 2 , "1", "" ),   # If the ratio count is >= 2, return the replicate number
                   ifelse( DT$Ratio_count_NCC_2 >= 2 , "2", "" ),
                   ifelse( DT$Ratio_count_NCC_3 >= 2 , "3", "" ))
DT[, Replicate_2rc_NCC := NCC[, paste(V1, V2, V3, sep = "") ]]     # Paste replicate numbers into one string and append to the big table

# For Rosco
Rosco <- data.table( ifelse( DT$Ratio_count_Rosco_1 >= 2 , "1", "" ),   # If the ratio count is >= 2, return the replicate number
                     ifelse( DT$Ratio_count_Rosco_2 >= 2 , "2", "" ),
                     ifelse( DT$Ratio_count_Rosco_3 >= 2 , "3", "" ))
DT[, Replicate_2rc_Rosco := Rosco[, paste(V1, V2, V3, sep = "") ]]     # Paste replicate numbers into one string and append to the big table

# For RoscoCPT
RoscoCPT <- data.table( ifelse( DT$Ratio_count_RoscoCPT_1 >= 2 , "1", "" ),   # If the ratio count is >= 2, return the replicate number
                        ifelse( DT$Ratio_count_RoscoCPT_2 >= 2 , "2", "" ),
                        ifelse( DT$Ratio_count_RoscoCPT_3 >= 2 , "3", "" ))
DT[, Replicate_2rc_RoscoCPT := RoscoCPT[, paste(V1, V2, V3, sep = "") ]]     # Paste replicate numbers into one string and append to the big table

# For TSA
TSA <- data.table( ifelse( DT$Ratio_count_TSA_1 >= 2 , "1", "" ),   # If the ratio count is >= 2, return the replicate number
                   ifelse( DT$Ratio_count_TSA_2 >= 2 , "2", "" ),
                   ifelse( DT$Ratio_count_TSA_3 >= 2 , "3", "" ))
DT[, Replicate_2rc_TSA := TSA[, paste(V1, V2, V3, sep = "") ]]     # Paste replicate numbers into one string and append to the big table


#### Combine the ATMi experiment ####

# Get mean and SEM over all six ATMi replicates at >= 2 ratio counts

# Get all six ratios with >= 2 ratio counts
R1 <- ifelse( DT$Ratio_count_ATMi_A_1 >= 2 , DT$log2_SILAC_ratio_ATMi_A_1 , NA  )
R2 <- ifelse( DT$Ratio_count_ATMi_A_2 >= 2 , DT$log2_SILAC_ratio_ATMi_A_2 , NA  )
R3 <- ifelse( DT$Ratio_count_ATMi_A_3 >= 2 , DT$log2_SILAC_ratio_ATMi_A_3 , NA  )

R4 <- ifelse( DT$Ratio_count_ATMi_B_1 >= 2 , DT$log2_SILAC_ratio_ATMi_B_1 , NA  )
R5 <- ifelse( DT$Ratio_count_ATMi_B_2 >= 2 , DT$log2_SILAC_ratio_ATMi_B_2 , NA  )
R6 <- ifelse( DT$Ratio_count_ATMi_B_3 >= 2 , DT$log2_SILAC_ratio_ATMi_B_3 , NA  )

# Combine into a new table
ATMi_comb <- data.table(R1, R2, R3, R4, R5, R6)

# How many ratios do we have for each protein?
ATMi_comb[, Number_of_ATMi_replicates_2rc := apply(.SD, 1, function(x){ sum(!is.na(x)) })]

# Get mean and SEM for proteins that were detected in at least 3 replicates 
# (ratios have already been median-normalised)
ATMi_ratio_cols <- c("R1", "R2", "R3", "R4", "R5", "R6")
ATMi_comb[, ATMi_combined_mean_ratio_2rc := apply(.SD, 1, mean, na.rm = TRUE), .SDcols = ATMi_ratio_cols]
ATMi_comb[, ATMi_combined_SEM_2rc        := apply(.SD, 1,   sd, na.rm = TRUE) / sqrt(Number_of_ATMi_replicates_2rc), .SDcols = ATMi_ratio_cols ]

# Put results back into DT
DT <- data.table( DT,
                  ATMi_comb[,.(ATMi_combined_mean_ratio_2rc,
                               ATMi_combined_SEM_2rc,
                               Number_of_ATMi_replicates_2rc)]  )


#### Add the sum of intensities per experiment ####

DT[, mean_intensity_ATMi_A := apply(.SD, 1, mean), .SDcols = grep("Intensity_ATMi_A", names(DT), value = TRUE) ]
DT[, mean_intensity_ATMi_B := apply(.SD, 1, mean), .SDcols = grep("Intensity_ATMi_B", names(DT), value = TRUE) ]
DT[, mean_intensity_CHX    := apply(.SD, 1, mean), .SDcols = grep("Intensity_CHX",    names(DT), value = TRUE) ]
DT[, mean_intensity_CPT    := apply(.SD, 1, mean), .SDcols = grep("Intensity_CPT",    names(DT), value = TRUE) ]
DT[, mean_intensity_HU     := apply(.SD, 1, mean), .SDcols = grep("Intensity_HU",     names(DT), value = TRUE) ]
DT[, mean_intensity_NCC    := apply(.SD, 1, mean), .SDcols = grep("Intensity_NCC",    names(DT), value = TRUE) ]
DT[, mean_intensity_Rosco  := apply(.SD, 1, mean), .SDcols = grep("Intensity_Rosco_",  names(DT), value = TRUE) ]
DT[, mean_intensity_RoscoCPT := apply(.SD, 1, mean), .SDcols = grep("Intensity_RoscoCPT", names(DT), value = TRUE) ]
DT[, mean_intensity_TSA      := apply(.SD, 1, mean), .SDcols = grep("Intensity_TSA",      names(DT), value = TRUE) ]
DT[, mean_intensity_ATMi_AB  := apply(.SD, 1, mean), .SDcols = grep("Intensity_ATMi",     names(DT), value = TRUE) ]


#### Write out results ####

fwrite(x = DT, file = "NCC9c_results_Apr19.csv")









