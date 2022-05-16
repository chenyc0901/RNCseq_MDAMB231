library(tidyverse)
library(edgeR)

rm(list=ls())
# Load data
lib_mtx <- read.csv("rawdata/library_info.csv", row.names = 1)
bsj_mtx <- read.csv("rawdata/circRNA_bsj.csv", row.names = 1, check.names=FALSE)
info_mtx <- read.csv("rawdata/circRNA_info.csv",check.names=FALSE)
MCF7_mtx <- bsj_mtx %>% 
            dplyr::select("RNA-MCF7","RNC-MCF7") %>% 
            filter(`RNA-MCF7`>3&`RNC-MCF7`>3)

#set 0.01 for technical replicates
bcv <- 0.01

circ_DElist <- DGEList(counts=MCF7_mtx,group=1:2)
dgeTest <- exactTest(circ_DElist,dispersion=bcv^2)

#extract the p value correction data
dgeTest_BH <- topTags(dgeTest,n=nrow(dgeTest$table))

#Draw volcanoplot
volcanoData <- cbind(dgeTest_BH$table$logFC,-log10(dgeTest_BH$table$FDR))
colnames(volcanoData) <- c("logFC","negLogPval")
plot(volcanoData,pch=0)

dgeTest_BH_table <- dgeTest_BH$table %>% rownames_to_column("id") 
MCF7_mtx_table <- MCF7_mtx %>% rownames_to_column("id")
MCF7_RNC_result <- left_join(dgeTest_BH_table,MCF7_mtx_table,by="id") 

#=================================================================================================================================================
 
MDAMB231_mtx <- bsj_mtx %>% 
  dplyr::select("RNA-MDAMB231","RNC-MDAMB231") %>% 
  filter(`RNA-MDAMB231`>3&`RNC-MDAMB231`>3)

#set 0.01 for technical replicates
bcv <- 0.01

circ_DElist <- DGEList(counts=MDAMB231_mtx,group=1:2)
dgeTest <- exactTest(circ_DElist,dispersion=bcv^2)

#extract the p value correction data
dgeTest_BH <- topTags(dgeTest,n=nrow(dgeTest$table))

#Draw volcanoplot
volcanoData <- cbind(dgeTest_BH$table$logFC,-log10(dgeTest_BH$table$FDR))
colnames(volcanoData) <- c("logFC","negLogPval")
plot(volcanoData,pch=0)

dgeTest_BH_table <- dgeTest_BH$table %>% rownames_to_column("id") 
MDAMB231_mtx_table <- MDAMB231_mtx %>% rownames_to_column("id")
MDAMB231_RNC_result <- left_join(dgeTest_BH_table,MDAMB231_mtx_table,by="id")
#==================================================================================================================================================
M7vs231_mtx <- bsj_mtx %>% 
  dplyr::select("RNC-MCF7","RNC-MDAMB231") %>% 
  filter(`RNC-MCF7`>3&`RNC-MDAMB231`>3)

#set 0.01 for technical replicates
bcv <- 0.01

circ_DElist <- DGEList(counts=M7vs231_mtx,group=1:2)
dgeTest <- exactTest(circ_DElist,dispersion=bcv^2)

#extract the p value correction data
dgeTest_BH <- topTags(dgeTest,n=nrow(dgeTest$table))

#Draw volcanoplot
volcanoData <- cbind(dgeTest_BH$table$logFC,-log10(dgeTest_BH$table$FDR))
colnames(volcanoData) <- c("logFC","negLogPval")
plot(volcanoData,pch=0)

dgeTest_BH_table <- dgeTest_BH$table %>% rownames_to_column("id") 
M7vs231_mtx_table <- M7vs231_mtx %>% rownames_to_column("id")
M7vs231_RNC_result <- left_join(dgeTest_BH_table,M7vs231_mtx_table,by="id")
#==========================================================================================================================================================

M7vs231_RNA_mtx <- bsj_mtx %>% 
  dplyr::select("RNA-MCF7","RNA-MDAMB231") %>% 
  filter(`RNA-MCF7`>3&`RNA-MDAMB231`>3)

#set 0.01 for technical replicates
bcv <- 0.01

circ_DElist <- DGEList(counts=M7vs231_RNA_mtx,group=1:2)
dgeTest <- exactTest(circ_DElist,dispersion=bcv^2)

#extract the p value correction data
dgeTest_BH <- topTags(dgeTest,n=nrow(dgeTest$table))

#Draw volcanoplot
volcanoData <- cbind(dgeTest_BH$table$logFC,-log10(dgeTest_BH$table$FDR))
colnames(volcanoData) <- c("logFC","negLogPval")
plot(volcanoData,pch=0)

dgeTest_BH_table <- dgeTest_BH$table %>% rownames_to_column("id") 
M7vs231_RNA_mtx_table <- M7vs231_RNA_mtx %>% rownames_to_column("id")
M7vs231_RNA_result <- left_join(dgeTest_BH_table,M7vs231_RNA_mtx_table,by="id")

#merge MCF7 MDAMB231 Result
Merged_result <- full_join(MCF7_RNC_result,MDAMB231_RNC_result,by="id") %>% 
                 full_join(.,M7vs231_RNA_result,by="id") %>%
                 full_join(.,M7vs231_RNC_result,by="id") %>% 
                 left_join(.,info_mtx,by=c("id"="circ_id")) %>%arrange(logFC.y.y) 

write_csv(Merged_result,file="edgeR_MDAMB231MCF7_RNCseq.csv")                                    
