# Load required libraries
library(data.table); library(ggplot2); library(grid); library(gridExtra)

# Read in the NCC data
NCC <- fread("NCC9_results_Aug17.csv")
setkey(NCC, SimpleID)


#### Contingencies ####

# Simplify lab internal categories
NCC$Lab_internal_category <- gsub("DNA repair",      "DNA repair and replication", NCC$Lab_internal_category)
NCC$Lab_internal_category <- gsub("DNA replication", "DNA repair and replication", NCC$Lab_internal_category)
NCC$Lab_internal_category <- gsub("Splicing",            "Splicing and RNA processing", NCC$Lab_internal_category)
NCC$Lab_internal_category <- gsub("pre-mRNA processing", "Splicing and RNA processing", NCC$Lab_internal_category)
NCC$Lab_internal_category <- gsub("pre-rRNA processing", "Splicing and RNA processing", NCC$Lab_internal_category)
NCC$Lab_internal_category <- gsub("Chromatin remodelling & assembly",   "Other chromatin factors", NCC$Lab_internal_category)
NCC$Lab_internal_category <- gsub("Chromatin structure & organization", "Other chromatin factors", NCC$Lab_internal_category)
NCC$Lab_internal_category <- gsub("DNA modifying",                      "Other chromatin factors", NCC$Lab_internal_category)
NCC$Lab_internal_category <- gsub("Histone-modifying activity",         "Other chromatin factors", NCC$Lab_internal_category)
NCC$Lab_internal_category <- gsub("Histones & histone variants",        "Other chromatin factors", NCC$Lab_internal_category)
NCC$Lab_internal_category <- gsub("Uncharacterized", "Uncharacterised - not annotated", NCC$Lab_internal_category)
NCC$Lab_internal_category[ is.na(NCC$Lab_internal_category) ] <- "Uncharacterised - not annotated"

NCC$Lab_internal_category <- factor(NCC$Lab_internal_category,
                                    levels=c("Uncharacterised - not annotated", "No expected chromatin function",
                                             "Transcription factor", "Splicing and RNA processing", "Other chromatin factors",
                                             "DNA repair and replication"))


# Make a simple contingency table (for 3+ ratio counts)
counts <- NCC[, c( grep("Lab_internal_category", names(NCC)), grep("Ratio_count", names(NCC))), with=FALSE] # Get the necessary columns

proteins_identified_N <- counts[, lapply(.SD, function(x){sum(x>=3)} ), .SDcols = !"Lab_internal_category"]
proteins_identified_N <- data.frame( proteins_identified_N = as.integer(proteins_identified_N),
                                     Replicates = gsub("Ratio_count_", "", names(proteins_identified_N)))
proteins_identified_N$Experiment <- gsub("_[[:digit:]]", "", proteins_identified_N$Replicates)

proteins_identified_N$Replicates <- factor(proteins_identified_N$Replicates, 
                                           levels = c("NCC_1", "NCC_2", "NCC_3", "Rosco_1", "Rosco_2", "Rosco_3", 
                                                      "CPT_1", "CPT_2", "CPT_3", "RoscoCPT_1", "RoscoCPT_2", "RoscoCPT_3", 
                                                      "ATMi_A_1", "ATMi_A_2", "ATMi_A_3", "ATMi_B_1", "ATMi_B_2", "ATMi_B_3",  
                                                      "HU_1", "HU_2", "HU_3",
                                                      "CHX_1", "CHX_2", "CHX_3", "TSA_1", "TSA_2", "TSA_3"))  # Set plotting order

pCon1 <- ggplot(proteins_identified_N, aes(x=Replicates, y=proteins_identified_N, fill=Experiment))+
         geom_bar(stat="identity")+
         geom_text(aes(label=proteins_identified_N, colour=Experiment), vjust=-0.5, size=2.5)+
         ylab("Number of identified proteins")+
         scale_fill_brewer(palette = "Set1")+
         scale_colour_brewer(palette = "Set1")+
         theme(panel.background = element_blank(), axis.title.x = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.25),
               axis.text.x=element_text(angle=90, hjust=1, size=6), axis.text.y=element_text(size=6), axis.ticks.x = element_blank(), 
               axis.ticks.y = element_line(size=0.25), axis.title.y = element_text(size=7),
               legend.text = element_text(size=6), legend.title = element_text(size=7),
               axis.line.x = element_line(size=0.25), axis.line.y = element_line(size=0.25))


# Break the contingency table down by lab internal category
identified_per_cat <- counts[, lapply(.SD, function(x){sum(x>=3)}) , by = Lab_internal_category] 
identified_per_cat <- melt(identified_per_cat, id.vars = "Lab_internal_category")
identified_per_cat$variable <- gsub("Ratio_count_", "", identified_per_cat$variable)
                
identified_per_cat$variable <- factor(identified_per_cat$variable, 
                                      levels = c("NCC_1", "NCC_2", "NCC_3", "Rosco_1", "Rosco_2", "Rosco_3", 
                                                 "CPT_1", "CPT_2", "CPT_3", "RoscoCPT_1", "RoscoCPT_2", "RoscoCPT_3", 
                                                 "ATMi_A_1", "ATMi_A_2", "ATMi_A_3", "ATMi_B_1", "ATMi_B_2", "ATMi_B_3",  
                                                 "HU_1", "HU_2", "HU_3",
                                                 "CHX_1", "CHX_2", "CHX_3", "TSA_1", "TSA_2", "TSA_3"))  # Set plotting order

pCon2 <- ggplot(identified_per_cat, aes(x=variable, y=value, fill=Lab_internal_category))+
         geom_bar(stat="identity")+
         ylab("Number of identified proteins")+
         scale_fill_manual(values=c("grey70", "firebrick2", "slateblue2", "aquamarine3", "deepskyblue2", "navy"))+
         theme(panel.background = element_blank(), axis.title.x = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.25),
               axis.text.x=element_text(angle=90, hjust=1, size=6), axis.text.y=element_text(size=6), axis.ticks.x = element_blank(), 
               axis.ticks.y = element_line(size=0.25), axis.title.y = element_text(size=7),
               legend.text = element_text(size=6), legend.title = element_text(size=7),
               axis.line.x = element_line(size=0.25), axis.line.y = element_line(size=0.25))


# Align and save the plots
pCon1 <- ggplotGrob(pCon1)
pCon2 <- ggplotGrob(pCon2)
pCon  <- rbind(pCon1, pCon2, size = "last")
grid.newpage()
grid.draw(pCon)

ggsave("contingency_plots_NCC9.pdf",
       pCon, height = 15, width = 19, units = "cm")

rm( list = ls()[! ls() %in% c("NCC")] )   # Clear workspace


#### Relationship between NCC experiments: correlation plots ####

# Add manual annotation from Martina
ReRe <- fread("Replication_repair_factors.csv")
ReRe <- ReRe[ !duplicated(ReRe$SimpleID) ]    # Remove proteins annotated more than once
setkey(ReRe, SimpleID)
NCC <- merge(NCC, ReRe, all.x = TRUE)

# Set ratios to NA unless they are based on >= 3 ratio counts
NCC$log2_SILAC_ratio_ATMi_A_1 <- ifelse( NCC$Ratio_count_ATMi_A_1 >= 3 , NCC$log2_SILAC_ratio_ATMi_A_1 , NA  )
NCC$log2_SILAC_ratio_ATMi_A_2 <- ifelse( NCC$Ratio_count_ATMi_A_2 >= 3 , NCC$log2_SILAC_ratio_ATMi_A_2 , NA  )
NCC$log2_SILAC_ratio_ATMi_A_3 <- ifelse( NCC$Ratio_count_ATMi_A_3 >= 3 , NCC$log2_SILAC_ratio_ATMi_A_3 , NA  )

NCC$log2_SILAC_ratio_ATMi_B_1 <- ifelse( NCC$Ratio_count_ATMi_B_1 >= 3 , NCC$log2_SILAC_ratio_ATMi_B_1 , NA  )
NCC$log2_SILAC_ratio_ATMi_B_2 <- ifelse( NCC$Ratio_count_ATMi_B_2 >= 3 , NCC$log2_SILAC_ratio_ATMi_B_2 , NA  )
NCC$log2_SILAC_ratio_ATMi_B_3 <- ifelse( NCC$Ratio_count_ATMi_B_3 >= 3 , NCC$log2_SILAC_ratio_ATMi_B_3 , NA  )

NCC$log2_SILAC_ratio_CHX_1 <- ifelse( NCC$Ratio_count_CHX_1 >= 3 , NCC$log2_SILAC_ratio_CHX_1 , NA  )
NCC$log2_SILAC_ratio_CHX_2 <- ifelse( NCC$Ratio_count_CHX_2 >= 3 , NCC$log2_SILAC_ratio_CHX_2 , NA  )
NCC$log2_SILAC_ratio_CHX_3 <- ifelse( NCC$Ratio_count_CHX_3 >= 3 , NCC$log2_SILAC_ratio_CHX_3 , NA  )

NCC$log2_SILAC_ratio_CPT_1 <- ifelse( NCC$Ratio_count_CPT_1 >= 3 , NCC$log2_SILAC_ratio_CPT_1 , NA  )
NCC$log2_SILAC_ratio_CPT_2 <- ifelse( NCC$Ratio_count_CPT_2 >= 3 , NCC$log2_SILAC_ratio_CPT_2 , NA  )
NCC$log2_SILAC_ratio_CPT_3 <- ifelse( NCC$Ratio_count_CPT_3 >= 3 , NCC$log2_SILAC_ratio_CPT_3 , NA  )

NCC$log2_SILAC_ratio_HU_1 <- ifelse( NCC$Ratio_count_HU_1 >= 3 , NCC$log2_SILAC_ratio_HU_1 , NA  )
NCC$log2_SILAC_ratio_HU_2 <- ifelse( NCC$Ratio_count_HU_2 >= 3 , NCC$log2_SILAC_ratio_HU_2 , NA  )
NCC$log2_SILAC_ratio_HU_3 <- ifelse( NCC$Ratio_count_HU_3 >= 3 , NCC$log2_SILAC_ratio_HU_3 , NA  )

NCC$log2_SILAC_ratio_NCC_1 <- ifelse( NCC$Ratio_count_NCC_1 >= 3 , NCC$log2_SILAC_ratio_NCC_1 , NA  )
NCC$log2_SILAC_ratio_NCC_2 <- ifelse( NCC$Ratio_count_NCC_2 >= 3 , NCC$log2_SILAC_ratio_NCC_2 , NA  )
NCC$log2_SILAC_ratio_NCC_3 <- ifelse( NCC$Ratio_count_NCC_3 >= 3 , NCC$log2_SILAC_ratio_NCC_3 , NA  )

NCC$log2_SILAC_ratio_Rosco_1 <- ifelse( NCC$Ratio_count_Rosco_1 >= 3 , NCC$log2_SILAC_ratio_Rosco_1 , NA  )
NCC$log2_SILAC_ratio_Rosco_2 <- ifelse( NCC$Ratio_count_Rosco_2 >= 3 , NCC$log2_SILAC_ratio_Rosco_2 , NA  )
NCC$log2_SILAC_ratio_Rosco_3 <- ifelse( NCC$Ratio_count_Rosco_3 >= 3 , NCC$log2_SILAC_ratio_Rosco_3 , NA  )

NCC$log2_SILAC_ratio_RoscoCPT_1 <- ifelse( NCC$Ratio_count_RoscoCPT_1 >= 3 , NCC$log2_SILAC_ratio_RoscoCPT_1 , NA  )
NCC$log2_SILAC_ratio_RoscoCPT_2 <- ifelse( NCC$Ratio_count_RoscoCPT_2 >= 3 , NCC$log2_SILAC_ratio_RoscoCPT_2 , NA  )
NCC$log2_SILAC_ratio_RoscoCPT_3 <- ifelse( NCC$Ratio_count_RoscoCPT_3 >= 3 , NCC$log2_SILAC_ratio_RoscoCPT_3 , NA  )

NCC$log2_SILAC_ratio_TSA_1 <- ifelse( NCC$Ratio_count_TSA_1 >= 3 , NCC$log2_SILAC_ratio_TSA_1 , NA  )
NCC$log2_SILAC_ratio_TSA_2 <- ifelse( NCC$Ratio_count_TSA_2 >= 3 , NCC$log2_SILAC_ratio_TSA_2 , NA  )
NCC$log2_SILAC_ratio_TSA_3 <- ifelse( NCC$Ratio_count_TSA_3 >= 3 , NCC$log2_SILAC_ratio_TSA_3 , NA  )


# Correlation plot for all proteins
ratios <- NCC[ , grep("log2_SILAC_ratio", names(NCC)), with = FALSE]    
names(ratios) <- gsub("log2_SILAC_ratio_", "", names(ratios))                       # Clean-up names

ratios <- ratios[, .( NCC_1, NCC_2, NCC_3, Rosco_1, Rosco_2, Rosco_3, CPT_1, CPT_2, CPT_3,
                      RoscoCPT_1, RoscoCPT_2, RoscoCPT_3, ATMi_A_1, ATMi_A_2, ATMi_A_3, ATMi_B_1, ATMi_B_2, ATMi_B_3,
                      HU_1, HU_2, HU_3, CHX_1, CHX_2, CHX_3, TSA_1, TSA_2, TSA_3)]  # Set plotting order

M <- cor(ratios, method = "pearson", use = "pairwise.complete.obs") # Get PCC
diag(M) <- NA                                                       # Remove self-comparisons
M <- melt(M, value.name = "PCC")                                    # Get shape ready for ggplot
M$Var2 <- factor(M$Var2, levels = rev(levels(M$Var2)) )             # Reverse factor levels for appropriate plotting order

pAll <- ggplot(M, aes(Var1, Var2, fill=PCC))+
        geom_tile()+
        geom_vline(xintercept = seq(0.5,27.5,3), size=0.25)+
        geom_hline(yintercept = seq(0.5,27.5,3), size=0.25)+
        ggtitle("Pearson's correlation (all proteins)")+
        scale_fill_gradient2( limits=c(-1,1), na.value = "grey80", low = "darkred", mid = "white", high = "mediumblue")+
        theme(panel.background = element_blank(), axis.title = element_blank(), axis.text.x=element_text(angle=90, hjust=1, size=6),
              axis.text.y=element_text(size=6), axis.ticks = element_blank(), legend.text = element_text(size=6),
              legend.title = element_text(size=7), plot.title = element_text(size=7, face="bold", hjust=0.5))


# Correlation plot restricted to DNA repair and replication factors
ratios <- NCC[ !is.na(Manual_annotation) , grep("log2_SILAC_ratio", names(NCC)), with = FALSE]    
names(ratios) <- gsub("log2_SILAC_ratio_", "", names(ratios))                       # Clean-up names

ratios <- ratios[, .( NCC_1, NCC_2, NCC_3, Rosco_1, Rosco_2, Rosco_3, CPT_1, CPT_2, CPT_3,
                      RoscoCPT_1, RoscoCPT_2, RoscoCPT_3, ATMi_A_1, ATMi_A_2, ATMi_A_3, ATMi_B_1, ATMi_B_2, ATMi_B_3,
                      HU_1, HU_2, HU_3, CHX_1, CHX_2, CHX_3, TSA_1, TSA_2, TSA_3)]  # Set plotting order

M <- cor(ratios, method = "pearson", use = "pairwise.complete.obs") # Get PCC
diag(M) <- NA                                                       # Remove self-comparisons
M <- melt(M, value.name = "PCC")                                    # Get shape ready for ggplot
M$Var2 <- factor(M$Var2, levels = rev(levels(M$Var2)) )             # Reverse factor levels for appropriate plotting order

pReRe <- ggplot(M, aes(Var1, Var2, fill=PCC))+
         geom_tile()+
         geom_vline(xintercept = seq(0.5,27.5,3), size=0.25)+
         geom_hline(yintercept = seq(0.5,27.5,3), size=0.25)+
         ggtitle("Pearson's correlation (DNA repair AND replication proteins)")+
         scale_fill_gradient2( limits=c(-1,1), na.value = "grey80", low = "darkred", mid = "white", high = "mediumblue")+
         theme(panel.background = element_blank(), axis.title = element_blank(), axis.text.x=element_text(angle=90, hjust=1, size=6),
               axis.text.y=element_text(size=6), axis.ticks = element_blank(), legend.text = element_text(size=6),
               legend.title = element_text(size=7), plot.title = element_text(size=7, face="bold", hjust=0.5))

pCor <- arrangeGrob(pAll, pReRe, nrow=2)   # Combine the plots
grid.newpage()
grid.draw(pCor)
ggsave("corr_plot_NCC9.pdf",
       pCor, height = 25, width = 15, units = "cm")

rm( list = ls()[! ls() %in% c("NCC")] )   # Clear workspace


#### Reproducibility: scatterplot between replicates ####

# Get ratio and annotation columns
ratios <- NCC[, c( grep("log2_SILAC_ratio", names(NCC), value = TRUE), "Manual_annotation"), with = FALSE]    
names(ratios) <- gsub("log2_SILAC_ratio_", "", names(ratios))          # Clean-up names

# Get all combinations between replicates
replicates <- names(ratios)[ grep("[[:digit:]]", names(ratios)) ]
replicates <- t( combn( replicates, 2))
replicates <- replicates[ gsub("_[[:digit:]]", "", replicates[,1]) == gsub("_[[:digit:]]", "", replicates[,2]) ,]

# Get r2 and plots for all replicate combinations
plot_list <- list()

for(i in 1:nrow(replicates)){
  exp_1 <- replicates[i,1]
  exp_2 <- replicates[i,2]
  temp  <- ratios[ , c(exp_1, exp_2, "Manual_annotation") , with=FALSE]
  colnames(temp)[1:2] <- c("x", "y")

PCC_all  <- cor( temp[,x], temp[,y], use = "pairwise.complete.obs" )
PCC_rere <- cor( temp[ !is.na(Manual_annotation) , x], temp[ !is.na(Manual_annotation) , y], use = "pairwise.complete.obs" )
r2_all  <- round( PCC_all*PCC_all   ,2)
r2_rere <- round( PCC_rere*PCC_rere ,2)

p <- ggplot(temp, aes(x,y))+
        geom_point(alpha=0.3, size=0.1, colour="grey50")+
        geom_smooth(method = "lm", se = FALSE, size=0.25, colour="grey50", fullrange=TRUE)+
        geom_point( data=temp[!is.na(Manual_annotation)] , alpha=0.8, size=0.5, colour="deepskyblue3")+
        geom_smooth(data=temp[!is.na(Manual_annotation)], method = "lm", se = FALSE, size=0.25, colour="deepskyblue3", fullrange=TRUE)+
        xlim(-2,2)+ylim(-2,2)+xlab(exp_1)+ ylab(exp_2)+
        annotate(geom = "text", x = -2, y = 2, label=paste("r2 =", r2_all),  vjust=1, hjust=0, size=2.5, colour="grey50")+
        annotate(geom = "text", x = 2, y = -2, label=paste("r2 =", r2_rere), vjust=0, hjust=1, size=2.5, colour="deepskyblue3")+
        theme(panel.background = element_blank(), axis.title.x = element_text(size=6, margin=margin(1,0,0,0)),
              axis.title.y = element_text(size=6, margin=margin(0,1,0,0)), axis.text = element_text(size=5), 
              axis.ticks = element_line(size=0.25), panel.border = element_rect(fill=NA, colour="black", size=0.25),
              axis.line.x = element_line(size=0.25), axis.line.y=element_line(size=0.25), plot.margin=unit(c(2,2,2,2), "mm"))

plot_list[[i]] <- p
}

# Assemble the plots into 1 figure (in appropriate order)
p <- arrangeGrob( plot_list[[16]], plot_list[[17]], plot_list[[18]],
                  plot_list[[19]], plot_list[[20]], plot_list[[21]],
                  plot_list[[10]], plot_list[[11]], plot_list[[12]],
                  plot_list[[22]], plot_list[[23]], plot_list[[24]],
                  plot_list[[1]], plot_list[[2]], plot_list[[3]],
                  plot_list[[4]], plot_list[[5]], plot_list[[6]],
                  plot_list[[13]], plot_list[[14]], plot_list[[15]],
                  plot_list[[7]], plot_list[[8]], plot_list[[9]], 
                  plot_list[[25]], plot_list[[26]], plot_list[[27]],
                  ncol=3)

ggsave("Replicates_NCC9.png", p,
       height = 30, width = 12, units = "cm")


rm( list = ls()[! ls() %in% c("NCC")] )   # Clear workspace


#### Histones ####

# Get all histones in dataset
Histones <- NCC[ Protein_names %like% "Histone H" ][ order(-NCC_Intensity) ]

# Simplify histone names and sort out histone types
Histones$SimpleHistoneName <- gsub(" type.+", "", Histones$Protein_names)
Histones$SimpleHistoneName <- gsub("H3.2", "H3", Histones$SimpleHistoneName)
Histones$SimpleHistoneName <- gsub(";Histone H1.0, N-terminally processed", "", Histones$SimpleHistoneName)
Histones$SimpleHistoneName <- gsub(";Histone ", " or\n ", Histones$SimpleHistoneName)
Histones$SimpleHistoneName <- gsub("Histone H3-like centromeric protein A", "CenpA", Histones$SimpleHistoneName)
Histones$SimpleHistoneName <- gsub("Histone ", "", Histones$SimpleHistoneName)

Histones[ SimpleHistoneName %in% c("H4", "H2B", "H3", "H2A"), HistoneType := "Core histones" ]
Histones[ SimpleHistoneName %like% "H1", HistoneType := "Linker histones" ]
Histones[ is.na(HistoneType) , HistoneType := "Histone variants" ]

# Drop H2A.J and H1t
Histones <- Histones[ !(SimpleHistoneName == "H2A.J" | SimpleHistoneName == "H1t") ]

# Turn into factors for plotting order
alphabetic_levels <- Histones[ order(HistoneType, SimpleHistoneName) ][ , SimpleHistoneName ]
Histones$SimpleHistoneName <- factor(Histones$SimpleHistoneName, levels = alphabetic_levels)

# Reshape data for plotting
ratios <- melt( Histones[, c( grep("mean_ratio", names(Histones), value = TRUE), "SimpleHistoneName", "HistoneType" ), with=FALSE ],
                id.vars = c("SimpleHistoneName", "HistoneType"), value.name = "mean_ratio")
sem <- melt( Histones[, c( grep("SEM", names(Histones), value = TRUE), "SimpleHistoneName", "HistoneType" ), with=FALSE ],
             id.vars = c("SimpleHistoneName", "HistoneType"), value.name = "SEM")
ratios$variable <- gsub("_mean_ratio", "", ratios$variable)
sem$variable <- gsub("_SEM", "", sem$variable)

histone_DT <- merge(ratios, sem, sort = FALSE)

# Histone barchart
pH1 <-  ggplot(histone_DT, aes(x=SimpleHistoneName, y=mean_ratio, fill=HistoneType))+
        geom_hline(yintercept = 0, colour="grey50", linetype="dashed")+
        geom_errorbar( aes(ymin = mean_ratio-SEM, ymax = mean_ratio+SEM , colour=HistoneType), size=0.25, width=0.5)+
        geom_bar(stat="identity")+
        facet_wrap(~variable, nrow = 4)+
        xlab("Histones")+
        ylab("mean SILAC ratio [log2]")+
        theme(panel.background = element_blank(), axis.title = element_text(size=7), axis.ticks.x = element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=6), axis.text.y=element_text(size=6),
              legend.text = element_text(size=6), legend.title = element_blank(),
              panel.border = element_rect(size=0.25, fill=NA, colour="black"), strip.background = element_rect(colour="black"))

ggsave("Histones_NCC9.pdf", pH1,
       height = 24, width = 19, units = "cm")


# Histone boxplot
histone_DT$Dummy_fill <- "Histones"   # Add a column to "dummy fill" the plot with same scale as the comparison plot later (pHB2)

pHB1 <- ggplot(histone_DT, aes(x=variable, y=mean_ratio, fill=Dummy_fill))+
        geom_hline(yintercept = 0, colour="grey50", linetype="dashed", size=0.25)+
        geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.2, size=0.25)+
        ylim(-1,1)+
        xlab("Experiment")+
        ylab("mean SILAC ratio [log2]")+
        scale_fill_brewer(palette = "Set1")+
        theme(panel.background = element_blank(), axis.title = element_text(size=7), axis.ticks.x = element_blank(),
              axis.text.x=element_text(size=6), axis.text.y=element_text(size=6),
              legend.text = element_text(size=6), legend.title = element_blank(),
              panel.border = element_rect(size=0.25, fill=NA, colour="black"), legend.background = element_blank(),
              legend.box.background = element_blank(), legend.key = element_blank())

# Histone boxplot compared to other protein groups
NCC[ Interphase_chromatin_probability > 0.6 , my_annotation := "Other chromatin"]
NCC[ Interphase_chromatin_probability < 0.3 , my_annotation := "Non-chromatin"]
NCC[ Majority_protein_IDs %in% Histones$Majority_protein_IDs , my_annotation := "Histones" ]

ratios <- NCC[, c( grep("mean_ratio", names(NCC), value = TRUE), "my_annotation"), with=FALSE ]
ratios <- ratios[ !is.na(my_annotation) ]

ratios <- melt(ratios, id.vars = "my_annotation", value.name = "mean_ratio")
ratios$variable <- gsub("_mean_ratio", "", ratios$variable)
ratios$my_annotation <- factor(ratios$my_annotation, levels=c("Histones", "Other chromatin", "Non-chromatin"))

pHB2 <- ggplot(ratios, aes(x=variable, y=mean_ratio, fill=my_annotation))+
        geom_hline(yintercept = 0, colour="grey50", linetype="dashed", size=0.25)+
        geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.2, size=0.25)+
        ylim(-1,1)+
        scale_fill_brewer(palette = "Set1")+
        xlab("Experiment")+
        ylab("mean SILAC ratio [log2]")+
        theme(panel.background = element_blank(), axis.title = element_text(size=7), axis.ticks.x = element_blank(),
              axis.text.x=element_text(size=6), axis.text.y=element_text(size=6),
              legend.text = element_text(size=6), legend.title = element_blank(),
              panel.border = element_rect(size=0.25, fill=NA, colour="black"), legend.background = element_blank(),
              legend.box.background = element_blank(), legend.key = element_blank())

# Align the histone boxplots
g1 <- ggplotGrob(pHB1)
g2 <- ggplotGrob(pHB2)
g  <- rbind(g1, g2, size="last")
grid.newpage()
grid.draw(g)

ggsave("Histone_boxplots_NCC9.pdf", g,
       height = 14, width = 19, units = "cm")


#### The new ATMi_B experiment ####

# Have a look at replication protein a and BRCA1 in the ATMi experiment

pA <- ggplot(NCC, aes(ATMi_B_mean_ratio, NCC_mean_ratio, label=Gene_names))+
       geom_point(alpha=0.3, size=0.2, colour="grey50")+
       geom_point(data = NCC[ Protein_names %like% "Replication protein A"] , col="red")+
       geom_point(data = NCC[ Gene_names == "BRCA1"] , col="blue")+
       xlim(-1,1)+
       ylim(-1,1)+
       ylab("NCC nascent - mature\nmean log2 SILAC ratio")+
       xlab("NCC - AMTi\nmean log2 SILAC ratio")+
       geom_text(data = NCC[ Protein_names %like% "Replication protein A"] , col="red", size=2, hjust=1.2, position=position_jitter(height = 0.05, width = 0))+
       geom_text(data = NCC[ Gene_names == "BRCA1"] , col="blue", size=2, hjust=-0.2)+
       geom_hline(yintercept = 0, size=0.25)+geom_vline(xintercept = 0, size=0.25)+
       theme(panel.background = element_blank(), axis.title.x = element_text(size=6, margin=margin(1,0,0,0)),
        axis.title.y = element_text(size=6, margin=margin(0,1,0,0)), axis.text = element_text(size=5), 
        axis.ticks = element_line(size=0.25), panel.border = element_rect(fill=NA, colour="black", size=0.25),
        axis.line.x = element_line(size=0.25), axis.line.y=element_line(size=0.25), plot.margin=unit(c(2,2,2,2), "mm"))

ggsave("NCC9_ATMi_vs_NM.pdf",
       pA, height = 12, width = 12, units = "cm")


# Reshape data for a barchart
candidates_ratios <- NCC[ Protein_names %like% "Replication protein A" | Gene_names == "BRCA1" , 
                          c("Gene_names", grep("_mean_ratio", names(NCC), value=T)), with=FALSE]
candidates_sem    <- NCC[ Protein_names %like% "Replication protein A" | Gene_names == "BRCA1" , 
                          c("Gene_names", grep("_SEM", names(NCC), value=T)), with=FALSE]

candidates_ratios <- melt(candidates_ratios, value.name = "mean_ratio", id.vars = "Gene_names")
candidates_ratios$variable <- gsub("_mean_ratio", "", candidates_ratios$variable)

candidates_sem <- melt(candidates_sem, value.name = "sem", id.vars = "Gene_names")
candidates_sem$variable <- gsub("_SEM", "", candidates_sem$variable)

candidates <- merge(candidates_ratios, candidates_sem)

pA2 <- ggplot(candidates, aes(x=variable, y=mean_ratio, fill=Gene_names))+
          geom_bar(stat="identity", position="dodge")+
          geom_errorbar( aes(ymin = mean_ratio-sem, ymax = mean_ratio+sem , colour=Gene_names), size=0.25, width=0.5, position=position_dodge(width=0.9))+
          geom_hline(yintercept = 0, size=0.25, colour="grey50", linetype="dashed")+
          scale_fill_brewer(palette = "Set1")+
          scale_colour_brewer(palette = "Set1")+
          ylab("mean SILAC ratio [log2]")+
          theme(panel.background = element_blank(), axis.title.y = element_text(size=7), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
            axis.text.x=element_text(angle=90, hjust=1, size=6), axis.text.y=element_text(size=6),
            legend.text = element_text(size=6), legend.title = element_blank(),
            panel.border = element_rect(size=0.25, fill=NA, colour="black"), strip.background = element_rect(colour="black"))
      
ggsave("NCC9_RPA_BRCA.pdf",
       pA2, height = 8, width = 19, units = "cm")


#### Other candidates Kyosuke wanted me to check ####

# (1) Repair
set1_ratios <- NCC[ Gene_names %in% c("FANCD2", "FANCI", "MDC1", "RAD18", "RNF168", "RNF169", "BLM") ,
                    c("Gene_names", grep("_mean_ratio", names(NCC), value=T)), with=FALSE]
set1_sem    <- NCC[ Gene_names %in% c("FANCD2", "FANCI", "MDC1", "RAD18", "RNF168", "RNF169", "BLM") ,
                    c("Gene_names", grep("_SEM", names(NCC), value=T)), with=FALSE]

set1_ratios <- melt(set1_ratios, value.name = "mean_ratio", id.vars = "Gene_names")
set1_ratios$variable <- gsub("_mean_ratio", "", set1_ratios$variable)

set1_sem <- melt(set1_sem, value.name = "sem", id.vars = "Gene_names")
set1_sem$variable <- gsub("_SEM", "", set1_sem$variable)

set1 <- merge(set1_ratios, set1_sem)

pS1 <- ggplot(set1, aes(x=variable, y=mean_ratio, fill=Gene_names))+
        geom_bar(stat="identity", position="dodge")+
        geom_errorbar( aes(ymin = mean_ratio-sem, ymax = mean_ratio+sem , colour=Gene_names), size=0.25, width=0.5, position=position_dodge(width=0.9))+
        geom_hline(yintercept = 0, size=0.25, colour="grey50", linetype="dashed")+
        scale_fill_brewer(palette = "Dark2")+
        scale_colour_brewer(palette = "Dark2")+
        ylab("mean SILAC ratio [log2]")+
        theme(panel.background = element_blank(), axis.title.y = element_text(size=7), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=6), axis.text.y=element_text(size=6),
              legend.text = element_text(size=6), legend.title = element_blank(),
              panel.border = element_rect(size=0.25, fill=NA, colour="black"), strip.background = element_rect(colour="black"))

ggsave("NCC9_repair.pdf",
       pS1, height = 8, width = 19, units = "cm")


# (2) PRC1
set2_ratios <- NCC[ Gene_names %in% c("BMI1", "RING1", "RNF2", "CBX2", "CBX4", "CBX8", "SCMH1", "PHC1", "PHC2", "PHC3") ,
                    c("Gene_names", grep("_mean_ratio", names(NCC), value=T)), with=FALSE]
set2_sem    <- NCC[ Gene_names %in% c("BMI1", "RING1", "RNF2", "CBX2", "CBX4", "CBX8", "SCMH1", "PHC1", "PHC2", "PHC3") ,
                    c("Gene_names", grep("_SEM", names(NCC), value=T)), with=FALSE]

set2_ratios <- melt(set2_ratios, value.name = "mean_ratio", id.vars = "Gene_names")
set2_ratios$variable <- gsub("_mean_ratio", "", set2_ratios$variable)

set2_sem <- melt(set2_sem, value.name = "sem", id.vars = "Gene_names")
set2_sem$variable <- gsub("_SEM", "", set2_sem$variable)

set2 <- merge(set2_ratios, set2_sem)

pS2 <- ggplot(set2, aes(x=variable, y=mean_ratio, fill=Gene_names))+
        geom_bar(stat="identity", position="dodge")+
        geom_errorbar( aes(ymin = mean_ratio-sem, ymax = mean_ratio+sem , colour=Gene_names), size=0.25, width=0.5, position=position_dodge(width=0.9))+
        geom_hline(yintercept = 0, size=0.25, colour="grey50", linetype="dashed")+
        scale_fill_brewer(palette = "Spectral")+
        scale_colour_brewer(palette = "Spectral")+
        ylab("mean SILAC ratio [log2]")+
        theme(panel.background = element_blank(), axis.title.y = element_text(size=7), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=6), axis.text.y=element_text(size=6),
              legend.text = element_text(size=6), legend.title = element_blank(),
              panel.border = element_rect(size=0.25, fill=NA, colour="black"), strip.background = element_rect(colour="black"))

ggsave("NCC9_PRC1.pdf",
       pS2, height = 8, width = 19, units = "cm")


# (3) H3K9me
set3_ratios <- NCC[ Gene_names %in% c("EHMT1", "EHMT2", "ZNF644") ,
                    c("Gene_names", grep("_mean_ratio", names(NCC), value=T)), with=FALSE]
set3_sem    <- NCC[ Gene_names %in% c("EHMT1", "EHMT2", "ZNF644") ,
                    c("Gene_names", grep("_SEM", names(NCC), value=T)), with=FALSE]

set3_ratios <- melt(set3_ratios, value.name = "mean_ratio", id.vars = "Gene_names")
set3_ratios$variable <- gsub("_mean_ratio", "", set3_ratios$variable)

set3_sem <- melt(set3_sem, value.name = "sem", id.vars = "Gene_names")
set3_sem$variable <- gsub("_SEM", "", set3_sem$variable)

set3 <- merge(set3_ratios, set3_sem)

pS3 <- ggplot(set3, aes(x=variable, y=mean_ratio, fill=Gene_names))+
        geom_bar(stat="identity", position="dodge")+
        geom_errorbar( aes(ymin = mean_ratio-sem, ymax = mean_ratio+sem , colour=Gene_names), size=0.25, width=0.5, position=position_dodge(width=0.9))+
        geom_hline(yintercept = 0, size=0.25, colour="grey50", linetype="dashed")+
        scale_fill_brewer(palette = "Set2")+
        scale_colour_brewer(palette = "Set2")+
        ylab("mean SILAC ratio [log2]")+
        theme(panel.background = element_blank(), axis.title.y = element_text(size=7), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=6), axis.text.y=element_text(size=6),
              legend.text = element_text(size=6), legend.title = element_blank(),
              panel.border = element_rect(size=0.25, fill=NA, colour="black"), strip.background = element_rect(colour="black"))

ggsave("NCC9_H3K9me.pdf",
       pS3, height = 8, width = 19, units = "cm")



