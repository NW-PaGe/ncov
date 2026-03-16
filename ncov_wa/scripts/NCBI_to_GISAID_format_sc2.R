#Taking NCBI WA data for formatting it for Nextstrain ingest


library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(tidyr)
library(stringr)
library(Biostrings)
#library(writexl)



print(.libPaths())
sessionInfo()
#getwd()

#Pass arguments from command line into this R script
args <- commandArgs(trailingOnly = TRUE)

input_metadata <- args[1]
input_sequences<- args[2]
output_metadata <- args[3]
output_sequences <- args[4]


##################################
# Processsing metadata file
##################################

#read in dataset

wa_metadata_raw <- read.delim(input_metadata)

clean_wa_metadata <- wa_metadata_raw %>%
        mutate(
          # Trim whitespace and standardize format
          Release_Date = str_trim(Isolate.Collection.date),
          date = str_trim(Isolate.Collection.date),
          date = str_replace_all(Isolate.Collection.date, "/", "-"),
          
          # Format Collection_Date to YYYY-MM-DD with "XX" placeholders
          date = case_when(
            str_detect(date, "^\\d{4}$") ~ paste0(date, "-XX-XX"),  # "2024" → "2024-XX-XX"
            str_detect(date, "^\\d{4}-\\d{2}$") ~ paste0(date, "-XX"),  # "2024-05" → "2024-05-XX"
            str_detect(date, "^\\d{4}-\\d{2}-\\d{2}$") ~ date,  # Keep full dates unchanged
            TRUE ~ NA_character_ ), # Keep other unknown formats as NA
          
          # Format Release_Date the same way
          Release_Date = case_when(
            str_detect(Release_Date, "^\\d{4}$") ~ paste0(Release_Date, "-XX-XX"),
            str_detect(Release_Date, "^\\d{4}-\\d{2}$") ~ paste0(Release_Date, "-XX"),
            str_detect(Release_Date, "^\\d{4}-\\d{2}-\\d{2}$") ~ Release_Date,
            TRUE ~ NA_character_),
          
          
          #Format Isolate.Lineage to strain
          strain=str_remove(Isolate.Lineage, pattern = "SARS-CoV-2/humans?/"),
          
        )%>%

    # Rename variables
        dplyr::rename(
          genbank_accession = Accession, 
          #date = Isolate.Collection.date, #this is done above
          host = Host.Name,
          region=Geographic.Region,
          sra_accession=SRA.Accessions,
          authors=Submitter.Names,
          pangolin_lineage=Virus.Pangolin.Classification
          ) %>%
    # Create Location variable
        mutate(Location = case_when(
            Geographic.Location %in% c("USA: Washington", "USA: WA") ~ "North America / United States / Washington",
            TRUE ~ Geographic.Location
          )) %>%
    # Create virus variable
        mutate(virus = "ncov" ) %>%
    # Extract Isolate_Name from GenBank_Title and clean up trailing subtype info
    #    mutate(Isolate_Name = str_extract(Isolate.Lineage, "\\(([^)]+)\\)")) %>%  
    #   mutate(Isolate_Name = str_remove_all(Isolate_Name, "[()]")) %>%        
    #    mutate(Isolate_Name = str_remove(Isolate_Name, "H[0-9]+N[0-9]+$")) %>%  
    # Format Geographic.Location to match GISAID format
      mutate(Geographic.Location_nospace = gsub("\\s*(:)\\s*", "\\1", Geographic.Location),
           location = "?") %>%
      mutate(location = gsub(".*?,\\s*", "", Geographic.Location)) %>%
      separate(Geographic.Location_nospace, into = c("country", "division"), sep = ":", fill = "right") %>%
      mutate(division = ifelse(is.na(division), "?", division)) %>%
    # Hard code division variable
    mutate(division = "Washington") 

#Write out metadata file subset to Washington sequences

write_tsv(clean_wa_metadata, output_metadata)

##################################
# Processsing fasta file
##################################
  
  
#  read FASTA files 

fasta <- readDNAStringSet(input_sequences, format = "fasta")


# Cleaning fasta names
#names(fasta)

name_list <- as.list(names(fasta))
fasta_names <- data.frame(name = unlist(name_list), stringsAsFactors = FALSE)

#names(fasta[,1])
names(fasta) <- sub(".*(USA/[^ ]+).*", "\\1", names(fasta) )


# Does metadata ID list match fasta sequence list?
meta_match <-clean_wa_metadata %>%
  filter(strain %in% names(fasta))

#Does fasta sequences match what is in metadata?
seq_match <- fasta[names(fasta) %in% meta_match$strain]


#Write to ncov directory
writeXStringSet(seq_match, filepath = output_sequences)

