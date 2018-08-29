setwd("~/Teagasc/all data rnaseq")

#loading libraries
if (FALSE){
  source("http://bioconductor.org/biocLite.R")
  biocLite("genefilter")
  biocLite("Rgraphviz")
  biocLite("ALL")
  biocLite("topGO")
  biocLite("multtest")
  biocLite("xtable")
  biocLite("globaltest")
  biocLite("xtable")
  biocLite("hgu133a.db")
  biocLite("hgu95av2.db")
  biocLite("ALL")
  biocLite("mvtnorm")
  biocLite("")
  biocLite("")
}

#biocLite("topGO")
if(TRUE){
  library('AnnotationDbi')
  library('graph')
  library("topGO")
  library("ALL")
  library("genefilter")
  library("Rgraphviz")
  library("ALL")
  library("topGO")
  library("multtest")
  library("xtable")
  library("globaltest")
  library("xtable")
  library("hgu133a.db")
  library("hgu95av2.db")
  library("mvtnorm")
  library("reshape2")
  library("reshape")
  library("tidyr")
  library("readxl")
  library("writexl")
  library("data.table")
  #library(stringr)
}

#Loading tables
  master<-as.data.table(read_xlsx("MASTER.xlsx"))
  filedata<-as.data.table(read.table("sleuth_wt_tableE_ABERGAIN")) 


#Renaming tables
  names(filedata)[1]<-"EVM_transcript"

#infos
  # Biological Process (BP):  
  # Cellular Component (CC):  	
  # Molecular Function (MF):  	
  

#Left outer join on EVM_transcript 
  filedata_merged<-master[filedata, on="EVM_transcript"]
  rm(filedata)
  
  
#### TopGO Analysis per see ####
  
#### Preparing data sets ####
#Preparing mapping GO file ggeneid2go.map  ####
  mapping<-master[,.("gene"=gene_id,"Go"=go)][!is.na(Go),][Go!="NA"]
  mapping2<-aggregate(Go ~ gene, mapping, FUN = toString) 
  
      #creating text file with all Go annotations for the genes in the "master" file
  write.table(mapping2, "ggeneid2go.map", sep="\t", row.names = FALSE,col.names = FALSE, quote = FALSE)
  rm(mapping)
  rm(mapping2)
    #importing mapping file for topGo
  geneID2GO <- readMappings(file = "ggeneid2go.map")
  #str(head(geneID2GO))

#making geneList for up and down regulated genes
#remove NA from the gene_id and pval in the data frame
#then remove duplicate
  geneL_up<-filedata_merged[(!is.na(gene_id)&!is.na(pval))&b>=0,][,.(gene_id,pval)]
  geneL_up<-geneL_up[!duplicated(geneL_up),]
  geneL_dn<-filedata_merged[(!is.na(gene_id)&!is.na(pval))&b<=0,][,.(gene_id,pval)]
  geneL_dn<-geneL_dn[!duplicated(geneL_dn),]


#making a Named num for up regulated genes, usable by topGo
  geneL_up_v<-as.vector(geneL_up$pval)
  names(geneL_up_v) <- geneL_up$gene_id
  rm(geneL_up)
  #str(geneL_up_v)

#making a Named num for dn regulated genes
  geneL_dn_v<-as.vector(geneL_dn$pval)
  names(geneL_dn_v) <- geneL_dn$gene_id
  rm(geneL_dn)


#preparing topgene list
#cutoff function

  topDiffGenes <- function(allScore) {
                       return(allScore < 0.01)
                        }

#number of selected genes
  xu <- topDiffGenes(geneL_up_v)
  xd <- topDiffGenes(geneL_dn_v)
  sum(xu)
  sum(xd)
  rm(xu)
  rm(xd)

#making the topGo object
  GOdata_up_BP <- new("topGOdata", 
              ontology = "BP", 
              allGenes = geneL_up_v,
              geneSel = topDiffGenes,
              annot = annFUN.gene2GO,
              nodeSize = 5,
              gene2GO = geneID2GO)

  GOdata_dn_BP <- new("topGOdata", 
                 ontology = "BP", 
                 allGenes = geneL_dn_v,
                 geneSel = topDiffGenes,
                 nodeSize = 5,
                 annot = annFUN.gene2GO,
                 gene2GO = geneID2GO)
  GOdata_up_BP
  GOdata_dn_BP

#### UP Regulated Genes Analysis #####
#all genes
  #a <- genes(GOdata_up_BP)
  #head(a)
  numGenes(GOdata_up_BP)
 
 #significan genes:
  #sg <- sigGenes(GOdata_up_BP)
  #str(sg)
  numSigGenes(GOdata_up_BP)
 
 #which GO terms are available for analysis and to obtain all the genes annotated to a subset of these
# GO terms
  #
  test.stat <- new("elimscore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
  resultELim_up_BP <- getSigGroups(GOdata_up_BP, test.stat)
  
  test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
  resultWeight_up_BP <- getSigGroups(GOdata_up_BP, test.stat)
  
  

   #based on count
    test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
    resultFisher_up_BP <- getSigGroups(GOdata_up_BP, test.stat)
   #based on score
    test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
    resultKS_up_BP <- getSigGroups(GOdata_up_BP, test.stat)
    
    test.stat <- new("weight01Score", testStatistic = GOKSTest, name = "KS tests")
    resultKSw01_up_BP <- getSigGroups(GOdata_up_BP, test.stat)
    
    test.stat <- new("weight01Score", testStatistic = GOKSTest, name = "GOglobalTest")
    resultglob_up_BP <- getSigGroups(GOdata_up_BP, test.stat)
    
   # whichAlgorithms()
    #whichTests()
    
    #head(resultKSw01_up_BP)
    #View(score(resultKSw01_up_BP ))
   # ?testStatistic()
    
    #allRes_up_BPt <- GenTable(GOdata_up_BP, classic = resultFisher_up_BP, Elim=resultELim_up_BP,
   #                          tresultKS_up_BP,  KS = resultKS_up_BP, weight = resultWeight_up_BP, a=resultglob_up_BP,
   #                          orderBy = "weight", ranksOf = "KS", topNodes = 200, numChar=100)
    
    #showSigOfNodes(GOdata_up_BP, score(resultglob_up_BP), firstSigNodes = 10, useInfo = 'all',.NO.CHAR = 100)
    
   # View(allRes_up_BPt)
    #Tests:
        if (TRUE)   { tresultFish <- runTest(GOdata_up_BP, algorithm = "classic", statistic = "fisher")
            tresultKS.elim_up_BP <- runTest(GOdata_up_BP, algorithm = "elim", statistic = "ks")
            tresultKS_up_BP <- runTest(GOdata_up_BP, algorithm = "classic", statistic = "ks")
            weight01.t <- runTest(GOdata_up_BP, algorithm = "weight01", statistic = "t")
            tallRes_up_BP <- GenTable(GOdata_up_BP, classicFisher = tresultFish,
                                       classicKS = tresultKS_up_BP, elimKS = tresultKS.elim_up_BP,
                                      weight01=weight01.t,
                                       orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
            showSigOfNodes(GOdata_up_BP, score(tresultKS_up_BP), firstSigNodes = 10, useInfo = 'all',.NO.CHAR = 100)
            View(tallRes_up_BP)  
            score(weight01.t)}
    
  #TEST MEthods
  # resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  # weight01.fisher <- runTest(GOdata, statistic = "fisher")
  # weight01.t <- runTest(GOdata, algorithm = "weight01", statistic = "t")
  # elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
   #head(score(resultWeight))
   #score(resultFisher)
    
   pvalFis_up_BP <- score(resultFisher_up_BP)
   #head(resultFisher)
   hist(pvalFis_up_BP, 50, xlab = "p-values")
   pvalWeight_up_BP <- score(resultWeight_up_BP, whichGO = names(pvalFis_up_BP))
   #head(pvalWeight)
   cor(pvalFis_up_BP, pvalWeight_up_BP)
   
   
   
   geneData(resultWeight_up_BP)
   allRes_up_BP <- GenTable(GOdata_up_BP, classic = resultFisher_up_BP, Elim=resultELim_up_BP,
                            tresultKS_up_BP,  KS = resultKS_up_BP, weight = resultWeight_up_BP,
                      orderBy = "weight", ranksOf = "KS", topNodes = 200, numChar=100)
   View(allRes_up_BP)
   #View(allRes)
   #Reformating the following columns as numeric
     allRes_up_BP$KS<-as.numeric(allRes_up_BP$KS)
     allRes_up_BP$classic<-as.numeric(allRes_up_BP$classic)
     allRes_up_BP$weight<-as.numeric(allRes_up_BP$weight)
     
    #environment(GOdata_up_BP@geneSelectionFun)[["filedata_merged"]][["description"]]
     
   #Export in xl file
   write_xlsx(allRes_up_BP,"analysis_up_BP.xlsx")

   #Making go2gene files
  # resultsgene2go<-data.table(go_id=character(), gene_list=list())
   
   #resultsgene2go<-rbind(resultsgene2go,c(goID,list(ann.genes)))
   #ann.genes <- genesInTerm(GOdata_up_BP, goID) ## get the annotations
   
   
   #View(resultsgene2go)
   
  # ann.genes <- genesInTerm(GOdata_up_BP, goID) ## get the annotations
  # head(ann.genes)
   
  #rm(resultsgene2go)  
    
    
  #7.3 Analysing individual GOs
  
   goID_up_BP <- allRes_up_BP[2, "GO.ID"]
   goID_up_BP <- "GO:0006412"
   print(showGroupDensity(GOdata_up_BP, goID_up_BP, ranks = TRUE))
   #str()
 
   #go_selected<-as.data.frame(allRes)
   #go_selected
   
   
  #Retriving genes

   
   #SHOW STRUCTURE
   par(cex = 0.3)
   showSigOfNodes(GOdata_up_BP, score(resultFisher_up_BP), firstSigNodes = 10, useInfo = 'all',.NO.CHAR = 100)
   showSigOfNodes(GOdata_up_BP, score(resultWeight_up_BP), firstSigNodes = 10, useInfo = 'def',.NO.CHAR = 100)
   
   showSigOfNodes(GOdata_up_BP, score(resultKS_up_BP), firstSigNodes = 10, useInfo = 'all',.NO.CHAR = 100)
   
   #printGraph(GOdata_up_BP, resultFisher_up_BP, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE,.NO.CHAR = 100)
   #printGraph(GOdata_up_BP, resultWeight_up_BP, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "def", pdfSW = TRUE,.NO.CHAR = 100)
  
  
  
  ### Retriving Term etc...
   graph(GOdata_up_BP)
   #show go term used
   ug <- usedGO(GOdata_up_BP)
   head(ug)
   goID_up_MF <- "GO:0006412"
   #We further select 10 random GO terms, count the number of annotated genes and obtain their annotation.
   sel.terms_up_BP <- "GO:0006412"# sample(usedGO(GOdata_up_BP), 10)
   
   num.ann.genes_up_BP <- countGenesInTerm(GOdata_up_BP, sel.terms_up_BP) ## the number of annotated genes
   
   show(ann.genes_up_BP)
   num.ann.genes_up_BP
   
   ann.genes_up_BP<- genesInTerm(GOdata_up_BP, sel.terms_up_BP) ## get the annotations
   head(ann.genes_up_BP)
   
   
   #ann.score <- scoresInTerm(GOdata_up_BP, sel.terms_up_BP)
   #head(ann.score_up_BP)
   #ann.score_up_BP <- scoresInTerm(GOdata_up_BP, sel.terms_up_BP, use.names = TRUE)
   #head(ann.score_up_BP)
   
   #termStat(GOdata_up_BP, sel.terms)
   
    
    
   #### DOWN Regulated Genes Analysis #####
   #all genes
   #ad <- genes(GOdata_dn_BP)
   #head(ad)
   numGenes(GOdata_dn_BP)
   
   #significan genes:
   #sgd <- sigGenes(GOdata_dn_BP)
   #str(sgd)
   numSigGenes(GOdata_dn_BP)
   
   #which GO terms are available for analysis and to obtain all the genes annotated to a subset of these
   # GO terms
   #
   test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
   resultWeightd_BP <- getSigGroups(GOdata_dn_BP, test.stat)
   
   #resultFis <- runTest(GOdata_dn_BP, algorithm = "classic", statistic = "fisher")
   #based on count
   test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
   resultFisherd_BP <- getSigGroups(GOdata_dn_BP, test.stat)
   #based on score
   test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
   resultKSd_BP <- getSigGroups(GOdata_dn_BP, test.stat)
   
   #TEST MEthods
   # resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
   # weight01.fisher <- runTest(GOdata, statistic = "fisher")
   # weight01.t <- runTest(GOdata, algorithm = "weight01", statistic = "t")
   # elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
   
   #head(score(resultWeightd_BP))
   #score(resultFisherd_BP)
   pvalFisd_BP <- score(resultFisherd_BP)
   #head(pvalFisd_BP)
   hist(pvalFisd_BP, 50, xlab = "p-values")
   pvalWeightd_BP <- score(resultWeightd_BP, whichGO = names(pvalFisd_BP))
   #head(pvalWeightd_BP)
   cor(pvalFisd_BP, pvalWeightd_BP)
   geneData(resultWeightd_BP)
   allResd_BP <- GenTable(GOdata_dn_BP, classic = resultFisherd_BP, KS = resultKSd_BP, weight = resultWeightd_BP,
                      orderBy = "weight", ranksOf = "KS", topNodes = 200, numChar=100)
   #View(allResd_BP)
   #Reformating the following columns as numeric
   allResd_BP$KS<-as.numeric(allResd_BP$KS)
   allResd_BP$classic<-as.numeric(allResd_BP$classic)
   allResd_BP$weight<-as.numeric(allResd_BP$weight)
   #Export in xl file
   write_xlsx(allResd_BP,"analysis_dn_BP.xlsx")
   
   #SHOW STRUCTURE
   showSigOfNodes(GOdata_dn_BP, score(resultFisherd_BP), firstSigNodes = 2, useInfo = 'all')
   showSigOfNodes(GOdata_dn_BP, score(resultWeightd_BP), firstSigNodes = 5, useInfo = 'def')
   
   showSigOfNodes(GOdata_dn_BP, score(resultKSd_BP), firstSigNodes = 5, useInfo = 'all',.NO.CHAR = 100)
   
   #printGraph(GOdata_dn_BP, resultFisherd_BP, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE,.NO.CHAR = 100)
   #printGraph(GOdata_dn_BP, resultWeightd_BP, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "def", pdfSW = TRUE,.NO.CHAR = 100)
   
   
   #gt <- printGenes(GOdata_dn_BP, whichTerms = goID, numChar = 40)
   
   #View(allRes)
   
   #environment(GOdata_dn_BP@geneSelectionFun)[["filedata_merged"]][["description"]]
   
 
   
   #Making go2gene files
   #resultsgene2god_BP<-data.table(go_id=character(),  gene_list=list())
   
   #resultsgene2go<-rbind(resultsgene2go,c(goID,list(ann.genes)))
   #ann.genes <- genesInTerm(GOdata_dn_BP, goID) ## get the annotations
   
   
   #View(resultsgene2go)
   
   #ann.genesd_BP <- genesInTerm(GOdata_dn_BP, goIDd_BP) ## get the annotations
   #head(ann.genesd_BP)
   
   #rm(resultsgene2god_BP)  
   
   
   #7.3 Analysing individual GOs
   
   goIDd_BP <- allResd_BP[1, "GO.ID"]
   print(showGroupDensity(GOdata_dn_BP, goIDd_BP, ranks = TRUE))
   #str(allRes)
   
   go_selected<-as.data.frame(allResd_BP)
   go_selected
   
   
   
   
   #### END #### 
   

#(sessionInfo())
  
   