#####################################
#####################################
###### trnL sequences analysis ######
#####################################
#####################################

# Before starting with the final filtering, we want to check the distribution of the number of times
#each sequence appears

# We filtered those sequences whose length is between 10 and 150 bp. 
library(ggplot2)
datos <- read.csv("C:/Users/CRISTINA/Downloads/distribucion_count.csv", sep = " ", header = F)
colnames(datos) <- c("Counts", "Sequences_number")
datos$Counts <- as.numeric(datos$Counts)
datos$Sequences_number <- as.numeric(datos$Sequences_number)

barplot(datos$Sequences_number, names.arg = datos$Counts,
        main = "Distribution of counts",
        xlab = "Counts", ylab = "Number of sequences occurring",
        col = "skyblue")

# Filtering those data that appears between 2 and 74 times in the whole dataset.
datos_filtrados <- subset(datos, Counts >= 2 & Counts <= 74)

barplot(
  height = datos_filtrados$Sequences_number,
  names.arg = datos_filtrados$Counts,
  main = "Distribution of counts (2-74 repeats)",
  xlab = "Counts",
  ylab = "Number of sequences",
  col = "skyblue",
  las = 1
)  

# % of occurences 
# First we calculate the number of total sequences we have in the dataset
total_secuencias <- sum(datos$Sequences_number)
# Add a percent column
datos$Percent <- round((datos$Sequences_number / total_secuencias) * 100, 2)
# Print the results
print(datos[, c("Counts", "Sequences_number", "Percent")])

datos_ordenados <- datos[order(-datos$Counts), ]

# Mean 
mean <- sum(datos$Counts * datos$Sequences_number) / sum(datos$Sequences_number)
print(paste("Man =", round(mean, 2)))
#Mean = 476.32

#Median
#Create a vector where each value of Counts is repeated according to its sequences_number.
vector_expandido <- rep(datos$Counts, datos$Sequences_number)
median <- median(vector_expandido)
print(paste("Median =", median))
#Median = 3


##############################################
##############################################
############ MOTUS & Occurrence ##############
##############################################
##############################################

library(dplyr)
library(stringr)
library(tibble)

# After the cleaning and filtering of the sequenceswith Obitools 4.4.0 we obtain 2 types of documents (csv). 
# The first document doesn´t have the occurrences. It contains the DNA sequence, their count,
# taxa, best database reference match, % of similarity...
#The second one has information about the occurrences of each sequence variant and the sample where we can find
# a specific sequence.The idea is to cmpare these 2 documents by merging the id column so we can have the 
# information of the taxa included in the occurrences file.

# Charge the first file
MOTUS_EMBL<- read.csv("C:/Users/CRISTINA/Downloads/MOTUS_trnL.csv", sep = ",")

#I create a new column called taxa taht contains the information of the column taxid but a little bit cleaned 
# without simbols...

EMBL <- MOTUS_EMBL %>%
  mutate(
    taxa = str_match(taxid, "\\[([^\\]]+)\\]")[,2]
  )

#The second csv contain the distribution of MOTU abundances across samples (correspond to individual PCR)

occurrences <- read.csv ("C:/Users/CRISTINA/Downloads/occurrency_0.1.csv")
#Traspose the table
nombres_columnas <- occurrences$id  # Guardar los CR_010P.1, CR_010P.10, etc.
occurrences <- occurrences[, -1, drop = FALSE] 
occur <- as.data.frame(t(occurrences))
colnames(occur) <- nombres_columnas
occur <- tibble::rownames_to_column(occur, var = "sample")


EMBL$id <- gsub(":", ".", EMBL$id)  
EMBL$id <- gsub("\\[", ".", EMBL$id)  
EMBL$id <- gsub("\\]", ".", EMBL$id)

#library(dplyr)
#Compare the id column with the sample column and creates a new column called taxa
occur <- occur %>%
  left_join(EMBL %>% select(id, taxa), 
            by = c("sample" = "id"))  # Join by coincidence

#Move the column taxa right after column 1 (sample)
occur <- occur %>%
  relocate(taxa, .after = 1)  

# There are manyyy sequences variants corresponding with multiple taxa. The main goal is to join the data by its
# taxa name to try to reduce the number of information. Each sequence variant has a taxa associated and a number
# that correspond with the occurences of that variant.
# If we join by the name of the taxa, we need to sum those values.
# So, next command is to join by taxa and sum the abundance column.
library(dplyr)

df_agrupado <- occur %>%
  group_by(taxa) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE),
            .groups = 'drop')

# We obtain 973 taxa

#Keep filtering to remove those data that don´t have information of taxa like: "cellular organism, 
# no rank or 50 kb inversion clade"

Filtered_EMBL_occurences <- df_agrupado %>%
  filter(
    !grepl("cellular organism|no rank|50 kb inversion clade", taxa, ignore.case = TRUE)
  )

write.csv(Filtered_EMBL_occurences, file = "EMBL_occurrences_filtered.csv", row.names = FALSE)

#PHYLOALPS database
#MOTUS_phyloAlps<- read.csv("C:/Users/CRISTINA/Downloads/MOTUS_PhyloAlps_10_COUNT.csv", sep = ",")



#GENERO UNA NUEVA COLUMNA LLAMADA TAXA PARA ELIMINAR LA "PORQUERIA" DE TAXID
EMBL <- Filtered_EMBL %>%
  mutate(
    taxa = str_match(taxid, "\\[([^\\]]+)\\]")[,2]
  )

Phylo <- Filtered_Alps %>%
  mutate(
    taxa = str_match(taxid, "\\[([^\\]]+)\\]")[,2]
  )

#Le cambio el nombre a las columnas "id" de ambas tablas
EMBL$id <- paste0("seq", 1:nrow(EMBL))

Phylo$id <- paste0("seq", 1:nrow(Phylo))







#16S.

anim <- read.csv("C:/Users/CRISTINA/Downloads/16S_motus.csv")

##GENERO UNA NUEVA COLUMNA LLAMADA TAXA PARA ELIMINAR LA "PORQUERIA" DE TAXID
animals <- anim %>%
  mutate(
    taxa = str_match(taxid, "\\[([^\\]]+)\\]")[,2]
  )


#occurrences
#genero una tabla trasposada a partir de la tabla de occurences

occ <- read.csv ("C:/Users/CRISTINA/Downloads/16S_occurrency.csv")

nombres_columnas <- occ$id  # Guardar los CR_010P.1, CR_010P.10, etc.
occ <- occ[, -1, drop = FALSE] 
occ <- as.data.frame(t(occ))
colnames(occ) <- nombres_columnas
occ_16S <- tibble::rownames_to_column(occ, var = "sample")


occ_16S$[,1] <- gsub(":", ".", occ_16S$[,1])  
occ_16S$id <- gsub("\\[", ".", occ_16S$id)  
occ_16S$id <- gsub("\\]", ".", occ_16S$id)

#library(dplyr)
#comparo id con los datos de sample y en base a eso agrego una nueva columna llamada taxa que contiene los taxones
occur <- %>%
  left_join(EMBL %>% select(id, taxa), 
            by = c("sample" = "id"))  # Une por coincidencias

#coloco la columna taxa en el segundo lugar despues del id de la muestra
occur <- occur %>%
  relocate(taxa.y, .after = 1)  # Mueve 'taxa' después de la 1ª columna


#taxa.y, manteniendo todas las demás columnas pero sumando los valores de abundancia cuando hay taxones repetidos.

# Agrupar por taxa.y y sumar las columnas de abundancia
df_agrupado <- occur %>%
  group_by(taxa.y) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE),
            .groups = 'drop')

#FILTRADO DE LAS SECUENCIAS DEL EMBL DATABASE
Filtered_EMBL <- df_agrupado %>%
  filter(
    !grepl("cellular organism|no rank|50 kb inversion clade", taxa.y, ignore.case = TRUE)
  )

write.csv(Filtered_EMBL, file = "EMBL_taxa_filtered.csv", row.names = FALSE)

#PHYLOALPS database
#MOTUS_phyloAlps<- read.csv("C:/Users/CRISTINA/Downloads/MOTUS_PhyloAlps_10_COUNT.csv", sep = ",")



#GENERO UNA NUEVA COLUMNA LLAMADA TAXA PARA ELIMINAR LA "PORQUERIA" DE TAXID
EMBL <- Filtered_EMBL %>%
  mutate(
    taxa = str_match(taxid, "\\[([^\\]]+)\\]")[,2]
  )

Phylo <- Filtered_Alps %>%
  mutate(
    taxa = str_match(taxid, "\\[([^\\]]+)\\]")[,2]
  )




#Le cambio el nombre a las columnas "id" de ambas tablas
EMBL$id <- paste0("seq", 1:nrow(EMBL))

Phylo$id <- paste0("seq", 1:nrow(Phylo))

