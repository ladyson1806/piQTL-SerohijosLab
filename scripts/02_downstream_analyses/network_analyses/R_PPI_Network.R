# Xavier Castellanos-Girouard
# 
# Date First Created: March 13 2023
# Date Last Updated: June 13 2024


# This script is dedicated to building a high-quality protein-protein
# interaction dataset.


#### Import libraries ####

library(tidyr)
library(dplyr)
library(org.Sc.sgd.db)
library(AnnotationDbi)

#### Import datasets ####

## Import IntAct dataset
# Full IntAct dataset (Release 243 - July 2022) was accessed on March 13 2023.
IntAct <-
  read.delim("/path/to/IntAct/intact.txt",
           sep = "\t",
           quote = "", 
           header = TRUE)


#View(ComplexPortal)

#### IntAct: Functions related to formatting ####

# This function allows for the extraction of loci names from the altID column
# of a dataframe containing IntAct data
getLocus <- 
  function(index, dataset, locus_col, alias_col){
    
    altID <- dataset[index, c(locus_col)]
    
    # Following regex function retrieves the 'ensemblfungi:' substrings and 
    # everything after (hopefully just the locus).
    locusID <-
      regmatches(altID, # Extract matched strings using positions
                 regexpr("ensemblfungi:.*", altID)) # Find positions of a match
    
    # Remove 'ensemblfungi:' substring from locusID
    locusID <- 
      sub(pattern = "ensemblfungi:",  
          replacement = "",
          locusID)
    
    # Remove any other substring from locusID
    locusID <- 
      sub(pattern = "\\|.*",  
          replacement = "",
          locusID)
    
    
    if (length(locusID) == 0){
      
      alias <- dataset[index, c(alias_col)]
      
      locusID <-
        regmatches(alias, 
                   regexpr("uniprotkb:......?.?.?.?.?\\(locus name\\)", alias))
      
      # Remove 'uniprotkb' substring
      locusID <- sub(pattern = "uniprotkb:", 
                    replacement = "",
                    locusID)
      
      # Remove '(gene name)' substring
      locusID <- sub(pattern = "(locus name)",
                    replacement = "",
                    locusID,
                    fixed = TRUE)
    }
    
    #print(paste(length(locusID), ":", locusID))
    
    # Return NA if no gene symbol was found, return symbol otherwise:
    if (length(locusID) == 0){return(NA)} else {return(locusID)}
  }

#### IntAct: Format and Filter IntAct data ####

## Only keep interactions where both interactors are from yeast.

# Boolean vector indicating whether interactorA is from yeast.
InteractorA_isYeast <-
  grepl(IntAct$Taxid.interactor.A, 
         pattern = "559292")

# Boolean vector indicating whether interactorB is from yeast.
InteractorB_isYeast <-
  grepl(IntAct$Taxid.interactor.B,
        pattern = "559292")

# Subset using boolean arrays
IntAct <- IntAct[(InteractorA_isYeast & InteractorB_isYeast),]


## Only keep interactions between proteins

# Boolean vector indicating whether interactorA is from yeast.
InteractorA_isProtein <-
  grepl(IntAct$Type.s..interactor.A, 
        pattern = "0326")

# Boolean vector indicating whether interactorB is from yeast.
InteractorB_isProtein <-
  grepl(IntAct$Type.s..interactor.B,
        pattern = "0326")

# Subset using boolean arrays
IntAct <- IntAct[(InteractorA_isProtein & InteractorB_isProtein),]


## Only keep interactions which were detected in yeast or in vitro

# list of unwanted hosts (i.e. anything not yeast or in vitro)
IntAct_host_remove <-
  c("taxid:83333(ecoli)|taxid:83333(\"Escherichia coli (strain K12)\")",
    "taxid:511693(ecobb)|taxid:511693(\"Escherichia coli (strain B / BL21)\")",
    "taxid:7111(trini)|taxid:7111(\"Trichoplusia ni (Cabbage looper)\")",
    "taxid:1772(mycsm)|taxid:1772(Mycobacterium smegmatis)",
    "taxid:9534(chlae-cos_7)|taxid:9534(Cercopithecus aethiops simian cells transformed with SV40)",
    "taxid:7108(spofr-sf_9)|taxid:7108(Spodoptera frugiperda insect cells)",
    "taxid:7108(spofr-sf_21)|taxid:7108(Spodoptera frugiperda insect cells)",
    "taxid:469008(ecobd)|taxid:469008(\"Escherichia coli (strain B / BL21-DE3)\")",
    "taxid:7111(trini-high_5)|taxid:7111(Trichoplusia ni cell line from eggs)",
    "taxid:562(ecolx)|taxid:562(Escherichia coli)",
    "taxid:9606(human-293t)|taxid:9606(Homo sapiens 293 cells transformed with SV40 large T antigen)",
    "taxid:7227(drome)|taxid:7227(\"Drosophila melanogaster (Fruit fly)\")",
    "taxid:7108(spofr)|taxid:7108(\"Spodoptera frugiperda (Fall armyworm)\")",
    "taxid:10090(mouse)|taxid:10090(Mus musculus)",
    "taxid:4934(lackl)|taxid:4934(Lachancea kluyveri)")

# Remove unwanted hosts
IntAct <- IntAct[!(IntAct$Host.organism.s. %in% IntAct_host_remove),]


## Remove interactions of type 'Association'

# Note: An Association is not rigorous enough for the purposes of the downstream
# analyses. (Also will save on computation time later on).

# Find interactions of type association
IntAct_isAssocation <- grepl(IntAct$Interaction.type.s., pattern = '0914')

# Remove interactions of type association
IntAct <- IntAct[!(IntAct_isAssocation),]

## Format columns

# Select only useful columns
IntAct <-
  IntAct %>%
  dplyr::select(X.ID.s..interactor.A, ID.s..interactor.B,
                Alt..ID.s..interactor.A, Alt..ID.s..interactor.B,
                Alias.es..interactor.A, Alias.es..interactor.B,
                Interaction.detection.method.s., Interaction.type.s.,
                Publication.1st.author.s., Publication.Identifier.s.,
                Confidence.value.s., Expansion.method.s.,
                Type.s..interactor.A, Type.s..interactor.B)


# Rename essential columns
colnames(IntAct)[1:6] <- 
  c("Uniprot_InteractorA", "Uniprot_InteractorB",
    "AltID_InteractorA", "AltID_InteractorB",
    "Alias_InteractorA", "Alias_InteractorB")


## Find gene name for each interactor

# For interactor A
IntAct$InteractorA <- 
  sapply(seq_along(IntAct$AltID_InteractorA),
        getLocus,
        dataset = IntAct,
        locus_col = "AltID_InteractorA",
        alias_col = "Alias_InteractorA")

# For interactor B
IntAct$InteractorB <- 
  sapply(seq_along(IntAct$AltID_InteractorB),
         getLocus,
         dataset = IntAct,
         locus_col = "AltID_InteractorB",
         alias_col = "Alias_InteractorB")


## Remove unmapped loci

# Find instances where value is not NA for InteractorA
InteractorA_notNA <- !is.na(IntAct$InteractorA)

# Find instances where value is not NA for InteractorB
InteractorB_notNA <- !is.na(IntAct$InteractorB)

# Only keep rows with no NA values
IntAct <- IntAct[(InteractorA_notNA & InteractorA_notNA),]


## Finish format with row index reset
row.names(IntAct) <- NULL

#### IntAct: Removal of duplicate interactions (First instance approach) ####

# This function determines whether an interaction is duplicated or not in a 
# dataset. Input is a character vector for genes (or their products), the second
# input is their interactors. Vectors must be of equal length, and is meant to 
# be organized such that the first element of vector A interacts with the 
# first element of vector B. Output is a boolean vector, containing a value for
# every interaction; TRUE if they are unique, FALSE otherwise (duplicated).
rm_dup_int <- function(InteractorA_vec, InteractorB_vec){
  
  top_i <- length(InteractorA_vec)
  
  is_duplicate_total <- c()
  # Find PPIs
  for (i in 1:(top_i-1) ){
    
    # Find all instances of an interactor in the vector
    # Only i+1 onwards needs to be check, eventually duplicated gene will return
    # FALSE; this is the one we keep.
    vecA_bool1 <- InteractorA_vec[i] == InteractorA_vec[(i+1):top_i] # For InteractorA
    vecB_bool1 <- InteractorB_vec[i] == InteractorB_vec[(i+1):top_i] # For InteractorB
    
    # If interactors have an instances at the same position, interaction
    # is present (TRUE element in vector).
    interaction_bool1 <- vecA_bool1 & vecB_bool1
    
    # Take into consideration Interaction AB and BA
    vecA_bool2 <- InteractorA_vec[i] == InteractorB_vec[(i+1):top_i] # For InteractorA
    vecB_bool2 <- InteractorB_vec[i] == InteractorA_vec[(i+1):top_i] # For InteractorB
    interaction_bool2 <- vecA_bool2 & vecB_bool2
    
    # Combine both
    interaction_bool <- interaction_bool1 | interaction_bool2
    
    # Check if the interaction is present in more than one place
    is_duplicate <- (sum(interaction_bool) > 0)
    
    is_duplicate_total <- c(is_duplicate_total, is_duplicate)
    }
  
  # Last one needs to be unique by necessity
  is_duplicate_total <- c(is_duplicate_total, FALSE)
  
  is_unique <- !is_duplicate_total
  
  return(is_unique)
}

IntAct_unique_bool <- rm_dup_int(IntAct$InteractorA, IntAct$InteractorB)
IntAct_unique1 <- IntAct[IntAct_unique_bool,]

IntAct_unique1 <- 
  IntAct_unique1[!duplicated(IntAct_unique1[, c("InteractorA", "InteractorB", "Publication.Identifier.s.")]),]

row.names(IntAct_unique1) <- NULL

IntAct_unique_tosend <- 
  IntAct_unique1[,c("Uniprot_InteractorA", "Uniprot_InteractorB",
                    "Confidence.value.s.", "InteractorA", "InteractorB")]


write.csv(IntAct_unique_tosend, "./data/IntAct_PPI_complete_polished.csv")




#### Export GeneName to ORF Table ####

### ORF To Gene Name


## Get bimap from org.Sc.sgd

x <- org.Sc.sgdGENENAME
# Get the gene names that are mapped to an ORF identifier
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])


## Extract ORFs from bimap

ORFs <- c()
Gene_names <- c()
for (ORF in names(xx)){
  
  if (length(xx[[ORF]] == 1)){
    Gene_names <- c(xx[[ORF]], Gene_names)
    ORFs <- c(ORF, ORFs)
  } else if (length(xx[[ORF]] >= 1)){
    for (i in xx[[ORF]]){
      Gene_names <- c(i, Gene_names)
      ORFs <- c(ORF, ORFs)
    }
  }
}

ORF_2_GeneName_DF <-
  data.frame(ORFs,
             Gene_names)

write.csv(ORF_2_GeneName_DF, "./data/Yeast_ORF2Gene.csv")
