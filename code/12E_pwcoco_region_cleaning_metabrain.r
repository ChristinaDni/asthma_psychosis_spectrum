library(dplyr)

setwd("/dir/revision_03_07_25/pwcoco/tissue/data/metabrain/raw") 

#Import the files
regions <- list.files(path = "/dir/revision_03_07_25/pwcoco/tissue/data/metabrain/raw", recursive = TRUE)

#Initiate a loop to start the cleaning
lapply(list(regions), function(k){
  
  for(k in regions) {
  
  tryCatch({
  
setwd("/dir/revision_03_07_25/pwcoco/tissue/data/metabrain/raw")

data<- vroom::vroom(k)
data<- as.data.frame(data)

#Create a filename column
#This is important to retain SNP information
data$filename<- basename(k)
data$filename<- gsub(".txt", "", data$filename)
print(head(data))

#Keep the file name- usefull for naming each region file
filename<- data$filename[1]
print(filename)

#Split each file by gene
map<- data %>% split(x=data, f=data$Gene_Ensembl)

  lapply(list(map), function(l){
  
  for(l in map) {
  
#Clean the number from the gene
l$clean_ensembl<- gsub("\\..*","",l$Gene_Ensembl)

#Keep the gene information
ensembl<- l$clean_ensembl[1]
print(ensembl)

#Keep the gene information
hgnc<- l$Gene_HGNC[1]
print(hgnc)

#Keep the columns necessary for pwcoco

region<- l%>%dplyr::select("SNP","A1","A2", "FRQ", "Beta", "SE", "P")
print(head(region))  

#Rename the column headers to match the asthma and all the rest
region<- region %>% dplyr::rename( 
     BETA = Beta
    ) 	

#Create a new name for the file
name<- paste0(filename, ":", ensembl, ":", hgnc, ".txt")

region$N<- ifelse(grepl("cerebellum", filename) == T, 492, 
			ifelse(grepl("cortex", filename) == T, 2683, 
			ifelse(grepl("hippocampus", filename) == T, 168, 
			ifelse(grepl("basalganglia", filename) == T, 208, 
			"NA"))))


#Save the file
setwd("/dir/revision_03_07_25/pwcoco/tissue/data/metabrain/clean")
write.table(region, name, sep="\t", quote=F, row.names=F, col.names=T)

}})
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}})
