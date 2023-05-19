ComputeBinaryMatrix <- function(gMCSAnalysis, gMCSList, genesHumanGEM, Cumulative = TRUE){
  
  List_gMCS_Essential_Matrix <- gMCSAnalysis$list.gMCS.essential.mat
  Binary_Matrix_Length_List <- vector("list", 1)
  
  Binary_Matrix <- matrix(0, nrow = length(genesHumanGEM),
                          ncol = ncol(gMCSAnalysis$mat.essential.gene),
                          dimnames = list(genesHumanGEM,
                                          colnames(gMCSAnalysis$mat.essential.gene)))
  
  for(Essential_Gene in names(List_gMCS_Essential_Matrix)){
      Binary_Matrix[Essential_Gene,] <- as.numeric(colSums(List_gMCS_Essential_Matrix[[Essential_Gene]]) > 0)

  }
  Binary_Matrix_Length_List<- Binary_Matrix
   
  return(Binary_Matrix_Length_List)
}