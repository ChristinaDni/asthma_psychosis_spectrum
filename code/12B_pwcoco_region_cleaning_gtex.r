library(dplyr)

setwd("/dir/revision_03_07_25/pwcoco/tissue/data/gtex/raw") 

#Import the files
regions <- list.files(path = "/dir/revision_03_07_25/pwcoco/tissue/data/gtex/raw", recursive = TRUE)

#Initiate a loop to start the cleaning
lapply(list(regions), function(k){
  
  for(k in regions) {
  
  tryCatch({
  
setwd("/dir/revision_03_07_25/pwcoco/tissue/data/gtex/raw")

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

#Split each file by probe
map<- data %>% split(x=data, f=data$Probe)

  lapply(list(map), function(l){
  
  for(l in map) {
  
#Clean the number from the probe
l$clean_probe<- gsub("\\..*","",l$Probe)

#Keep the probe information
probe<- l$clean_probe[1]
print(probe)

#Keep the gene information
gene<- l$Gene[1]
print(gene)

#Keep the columns necessary for pwcoco

region<- l%>%dplyr::select("SNP","A1","A2", "Freq", "b", "SE", "p")
print(head(region))  

#Rename the column headers to match the asthma
region<- region %>% dplyr::rename(
      FRQ = Freq,
     BETA = b,
      P = p
    ) 
	
#Add the sample size
region$N<- 444

#Create a new name for the file
name<- paste0(filename, ":", probe, ":", gene, ".txt")

#Save the file
setwd("/dir/revision_03_07_25/pwcoco/tissue/data/gtex/clean")
write.table(region, name, sep="\t", quote=F, row.names=F, col.names=T)

}})
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}})

###

