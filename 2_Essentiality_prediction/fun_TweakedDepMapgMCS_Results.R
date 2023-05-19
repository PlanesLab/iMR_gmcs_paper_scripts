TweakedDepMapgMCS_Results <- function(Binary_Matrix, 
                                      Effect,
                                      Table_HumanGEM_Genes,
                                      gMCS_order,
                                      Threshold_Info)
  
  
  try({
    result_list <- list()
    evalMethod = "All_Tasks"
    # browser()
    ExpMatBinary <- Effect < (-0.6)
    ExpMatBinary[is.na(Effect)] <- 0 
    tissues <- intersect(colnames(Binary_Matrix),
                         colnames(ExpMatBinary))
    
    # Initialize Variables #
    results = data.frame('cellLine' = tissues,
                         'TP' = NA, 'TN' = NA, 'FP' = NA, 'FN' = NA,
                         "check_sum_genes" = NA, 
                         'accuracy' = NA, 'sensitivity' = NA,
                         'specificity' = NA, 'F1' = NA,
                         'MCC' = NA, 'Penr' = NA, 'logPenr' = NA,
                         'PenrAdj' = NA, 'logPenrAdj' = NA)
    
    for (t in 1:length(tissues)){
      
      # Genes From the Model #
      modelGenes       <- unique(Table_HumanGEM_Genes$ENSEMBL) 
      modelGenesSYMBOL <- unique(Table_HumanGEM_Genes$SYMBOL) 
      # Essential Genes #
      modelEssential        <- rownames(Binary_Matrix)[Binary_Matrix[,tissues[t]]>0]
      modelEssentialENSEMBL <- unique(Table_HumanGEM_Genes$ENSEMBL[Table_HumanGEM_Genes$ENSEMBL %in% modelEssential])
      modelEssentialSYMBOL  <- unique(Table_HumanGEM_Genes$SYMBOL[Table_HumanGEM_Genes$ENSEMBL %in% modelEssential])
      # Non-Essential Genes #
      modelNonEssential       <- setdiff(modelGenes, modelEssential)
      modelNonEssentialSYMBOL <- unique(Table_HumanGEM_Genes$SYMBOL[Table_HumanGEM_Genes$ENSEMBL %in% modelNonEssential])
      modelNonEssentialENSEMBL <- unique(Table_HumanGEM_Genes$ENSEMBL[Table_HumanGEM_Genes$ENSEMBL %in% modelNonEssential])
      
      # expGenes has all the genes that are not NA from Achilles #
      expGenes <-  rownames(ExpMatBinary)[!is.na(ExpMatBinary[,tissues[t]])]
      # expEssential has all the genes that are considered essential #
      expEssential <- expGenes[ExpMatBinary[expGenes,tissues[t]]==1]
      # expNonEssential has all the genes that are NOT essential # 
      expNonEssential <- setdiff(expGenes,expEssential)
      
      
      results$TP[t] = length(intersect(modelEssentialENSEMBL, expEssential))        # True Positives  #
      results$TN[t] = length(intersect(modelNonEssentialENSEMBL, expNonEssential))  # True Negatives  #
      results$FP[t] = length(intersect(modelEssentialENSEMBL, expNonEssential))     # False Positives #
      results$FN[t] = length(intersect(modelNonEssentialENSEMBL, expEssential))     # False Negatives #   
      
      # Intersect Genes from HumanGEM and Genes from ACHILLES #
      pop <- intersect(modelGenesSYMBOL, expGenes)
      # Essential Genes from gMCST5 #
      sample <- intersect(modelEssentialSYMBOL, expGenes)
      # Essential Genes from ACHILLES #
      successes <- intersect(expEssential, modelGenesSYMBOL)
      
      results$Penr[t] <- phyper(length(intersect(successes,sample)) - 1, # Number of White Balls Drawn  #
                                length(intersect(successes,pop)),        # Number of White Balls in Urn #
                                length(setdiff(pop, successes)),         # Number of Black Balls in Urn #
                                length(sample),                          # Total Number of Balls Drawn  #  
                                lower.tail = FALSE);
      
    }
    
    
    results$check_sum_genes =  results$TP + results$TN + results$FP + results$FN
    
    
    
    # Calculate Metrics
    results$sensitivity = results$TP/(results$TP + results$FN)
    results$specificity = results$TN/(results$TN + results$FP)
    results$accuracy    = (results$TP + results$TN)/(results$TP + results$TN + results$FP + results$FN)
    results$PPV         = results$TP/(results$TP+results$FP)
    results$F1          = 2*results$TP/(2*results$TP + results$FP + results$FN)
    results$MCC         = ((results$TP*results$TN) - (results$FP*results$FN))/(sqrt((results$TP+results$FP)*(results$TP+results$FN))*sqrt((results$TN+results$FP)*(results$TN+results$FN)))  # Matthews correlation coefficient
    results$PenrAdj     = p.adjust(results$Penr,'BH')
    results$logPenr     = -log10(results$Penr)
    results$logPenrAdj  = -log10(results$PenrAdj)
    
    results$model <- "HumanGEM"
    results$evalMethod <- evalMethod
    
    results$thMethod <- Threshold_Info
    if(gMCS_order > 7){
      results$Order    <- paste0("HigherThan_Order_7")
    } else {
      results$Order    <- paste0("Order_", gMCS_order)
    }
    
    results$PenrAdj = p.adjust(results$Penr,'BH')
    
    results$logPenr = -log10(results$Penr)
    results$logPenrAdj = -log10(results$PenrAdj)
    
    result_list[[evalMethod]] <- results
    
    results <- do.call(rbind, result_list)
    return(results)})