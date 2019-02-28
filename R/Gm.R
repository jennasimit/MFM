#' Simulated genotype matrix as class SnpMatrix 
#'
#' Simulated data was generated using functions in MFMextra. In this example, there are two diseases that each have two
#' causal variants, of which one is shared beween diseases; both have causal variant rs61839660 in
#' group A; disease 1 has additional causal variant rs56382813 in group D; disease 2 has additional causal variant rs11594656 in
#' group C (SNP groups are given in the snpGroups object of this package). There are 3000 each of disease 1, disease 2, 
#' and controls. 
#'
#' This SnpMatrix has 9000 rows and 26 columns.
#' Rows are individuals and are labelled as control (1:3000), case1 (1:3000), or case2 (1:3000). 
#' Columns are SNPs and labelled by rsID
#'
"Gm"
