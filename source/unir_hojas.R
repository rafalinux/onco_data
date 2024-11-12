setwd("~/Python_Projects/ONCO_Genes")
hoja_1 <- read.csv("hoja_1.csv", header=TRUE, sep=",")
hoja_2 <- read.csv("hoja_2.csv", header=TRUE, sep=",")
hoja_2$tumor_id <- hoja_2$Sample_Name_Corregido
hojas_merge <- merge(hoja_1, hoja_2, by = "tumor_id", 
                  all.x = TRUE) 
library("tidyr")

# Remove rows with NA's using drop_na()
hojas_merge <- hojas_merge %>% drop_na()

hojas_merge$Sample_Name_Corregido <- NULL
hojas_merge$Pte. <- NULL

write.csv(hojas_merge, file = "data_onco.csv", row.names = TRUE)
