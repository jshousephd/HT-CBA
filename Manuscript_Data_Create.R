#############################################################################################################
##Author :          John House - BRC;NCSU
##Last Modified :   2/28/2017 - John House
##Purpose:          Subset Cardiomyocyte Data for manuscript to sample data set
##
#############################################################################################################

setwd("C:/Users/jshouse/Google_Drive/TempOSeq_Pipeline/ht_cba_manuscript")

### READ IN COUNT MATRIX AND ATTRIBUTE FILE (HASH FILE) ###
raw_count_data <- read.csv(file = "raw_counts.csv",header=T,stringsAsFactors = F, row.names = 1) %>% as.data.frame()
hash <- read.csv(file = "hashfile.csv",stringsAsFactors = FALSE,header=T,row.names=11)
head(hash)
manuscript_hash <- subset(hash,hash$Treatment %in% c("Isoproterenol","Propranolol","Nifedipine","dofetilide") |
                            (hash$Index.Set %in% c("B","C") & hash$Treatment %in% c("MEDIA","VEHICLE")))
manuscript_hash$HEADER <- paste0("X",rownames(manuscript_hash))
manuscript_counts <- raw_count_data[,manuscript_hash$HEADER]

write.csv(manuscript_hash,file="Input Files/hash.csv",row.names=FALSE)
write.csv(manuscript_counts,file="Input Files/raw_counts.csv")



      
    