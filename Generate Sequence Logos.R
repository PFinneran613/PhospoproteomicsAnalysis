plotLogo <- function(kinase){
  library("rjson")
  require(ggplot2)
  require(ggseqlogo)
  data(ggseqlogo_sample)
  
  AAs = c('I', 'Q', 'M', 'D', 'L', 'S', 'T', 'W', 'G','H', 'A','E', 'P', 'Y', 'V', 'K','C', 'R', 'N', 'F')
  colors = c('black', '#A97C50', 'black', '#800000', 'black', '#4A79A5', '#4A79A5', '#6F6F6F', '#1C5E3F','#142B4F', 'black','#800000', '#1C5E3F', '#6F6F6F', 'black', '#142B4F','#BDB76B', '#142B4F', '#A97C50', '#6F6F6F')
  Pos = c(-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7)
  fname1 <-sprintf("PhosphoProteomes/%s/%s_Scores.csv",kinase,kinase)
  data <- read.csv(file=fname1, head=FALSE,sep=",")
  Data <- as.matrix.data.frame(data,rownames.force=TRUE,rownames.value=AAs)
  rownames(Data) <- AAs
  colnames(Data) <- Pos
  
  cs1 = make_col_scheme(chars=AAs, cols=colors) 
  
  ggplot() + 
    annotate('rect', xmin = 7.5, xmax = 8.5, ymin = 0.0, ymax = Data[14,8], alpha = 1, fill='#CCCCCC') +
    geom_logo(Data, method='custom', seq_type='aa',col_scheme=cs1,rev_stack_order=TRUE) +
    theme_logo() + theme(axis.text.x = element_blank(),axis.text.y = element_blank())
  
  fname2 <- sprintf('PhosphoProteomes/%s/%s Logo.eps',kinase, kinase)
  ggsave(file=fname2, width=4.17, height=1.5,device="eps")
  fname2 <- sprintf('PhosphoProteomes/%s/%s Logo.png',kinase, kinase)
  ggsave(file=fname2, width=8.34, height=3,device="png")
}


here <- dirname(parent.frame(2)$ofile)
setwd(here)


kinases <- list('Abl','Src','Anc S1', 'Anc AS', 'Anc AST', 'Src PSP', 'Abl PSP', 'Src 10 Minute', 'Src 10 Minute B','Swiss Prot','MS Data Background')

for(kinase in kinases){plotLogo(kinase)}
  
