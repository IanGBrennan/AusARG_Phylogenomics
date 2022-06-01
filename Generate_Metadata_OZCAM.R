require(dplyr)
require(galah); 
require(googlesheets4)
require(stringr)
require(taxizedb)

# To query the ALA database through 'galah' you need to have an account with the ALA.org.au
## once you have an account you'll configure 'galah' using that email address: (or use mine)
ala_config(email="iangbrennan@gmail.com")

# Note on museum records:
## Preferentially use MUSEUM numbers over ABTC numbers
## registration numbers from the SAMA, WAM, MAGNT, CSIRO/ANWC begin RXXXXXX
## registration numbers from the AMS begin with R.XXXXXX
## registration numbers from the QM begin with JXXXXXX
## registration numbers from the MV begin with DXXXXXX or RXXXXXX
## registration numbers from the ABTC begin with ABTCXXXXXX

# remove authorization token
gs4_deauth()

# read in the file of samples to collect metadata for
# THIS REQUIRES COLUMNS: RegNo, Genus, Species (named exactly), Extraction
extracts <- read_sheet("https://docs.google.com/spreadsheets/d/1pswHDW2rq389WomTh5Crg59RPXDV9r6B160P5AV7q4Y/edit#gid=0") # pygos
#extracts <- read_sheet("https://docs.google.com/spreadsheets/d/15jPBNIGyvrhjUp0CWgw14ZqKVrHoQq1F8fPedD_-UA0/edit#gid=0") # Anilios
#extracts <- read_sheet("https://docs.google.com/spreadsheets/d/13NkCufHJMIbReiZeTXsXCsuRg0MFTrIXsur7pYNdKik/edit#gid=0") # diplos


## if your column RegNo includes the museum code (e.g. WAM R123456)
#extracts <- mutate(extracts, RegNo = sapply(extracts$RegNo, function(x) str_split(x, " ")[[1]][2]),
#                             institution = sapply(extracts$RegNo, function(x) str_split(x, " ")[[1]][1]))

# ABTC numbers are formatted differently, add 'ABTC' before the number
#extracts[which(extracts$institution == "ABTC"),"RegNo"] <- 
#  sapply(extracts[which(extracts$institution == "ABTC"),"RegNo"], function(x) paste0("ABTC",x))

# if you have samples you'd like to exclude from this activity
extracts <- filter(extracts, AusARG == "yes")

# pull the ala records for your broad group
records <- ala_occurrences(taxa = select_taxa("Pygopodidae"),
                           columns = select_columns("scientificName",
                                                    "decimalLatitude",
                                                    "decimalLongitude",
                                                    "institutionName",
                                                    "collectionName",
                                                    "institutionCode",
                                                    "catalogNumber",
                                                    "otherCatalogNumbers",
                                                    "phylum",
                                                    "class",
                                                    "order",
                                                    "family",
                                                    "genus",
                                                    "species",
                                                    "subspecies",
                                                    "common_name_and_lsid",
                                                    "country",
                                                    "stateProvince",
                                                    "verbatimLocality",
                                                    "habitat",
                                                    "eventDate",
                                                    "coordinateUncertaintyInMeters",
                                                    "preparations",
                                                    "sex",
                                                    "lifeStage",
                                                    "recordID",
                                                    #"taxonID",
                                                    "locality",
                                                    "collector",
                                                    "identifiedBy",
                                                    "typeStatus"))
# fix the locality field in the records
records <- mutate(records, locality = paste(locality, "|", verbatimLocality))

# fix the catalogNumber field (some ABTC records are ABTC1111.1)
ABTCs <- dplyr::filter(records, collectionName == "South Australian Museum Australian Biological Tissue Collection")
ABTCs$catalogNumber <- sapply(ABTCs$catalogNumber, function(x) strsplit(x, "[.]")[[1]][1])
records <- rbind(records, ABTCs) # combine the fixed ABTC with old data
records <- dplyr::distinct(records, catalogNumber, .keep_all=T) # remove duplicate records

# extract info for the samples that are in your extraction list
easy_matches <- filter(records, catalogNumber %in% extracts$RegNo); 
paste("there are", (nrow(extracts)-nrow(easy_matches)), "records to curate by hand")

# ALA has unreliable names, ours are probably better, so change accordingly
for (kk in 1:nrow(easy_matches)){
  print(kk)
  curr.match <- easy_matches[kk,]
  easy_matches[kk,"scientificName"] <- paste(extracts[which(extracts$RegNo == curr.match$catalogNumber),"Genus"],
                                             extracts[which(extracts$RegNo == curr.match$catalogNumber),"Species"])
  easy_matches[kk,"genus"] <- extracts[which(extracts$RegNo == curr.match$catalogNumber),"Genus"]
  easy_matches[kk,"species"] <- easy_matches[kk,"scientificName"]
  easy_matches[kk,"Extraction"] <- extracts[which(extracts$RegNo == curr.match$catalogNumber),"Extraction"]
}

# find the ones that slipped through the cracks (no ALA records)
#missed <- extracts[which(!extracts$number %in% records$catalogNumber),]
missed <- extracts[which(!extracts$RegNo %in% records$catalogNumber),]


# make empty data frame to add in the missed samples
missing <- data.frame(matrix(data = "", nrow = nrow(missed), ncol = ncol(easy_matches)))
colnames(missing) <- colnames(easy_matches)
missing <- mutate(missing, scientificName = paste(missed$Genus, missed$Species),
                           genus = missed$Genus, species = paste(missed$Genus, missed$Species),
                           #catalogNumber = missed$number, institutionCode = missed$institution)
                           catalogNumber = missed$RegNo, institutionCode = "?", Extraction = missed$Extraction)


# add the samples lacking ALA data to the bottom of your successful searches
all_matches <- rbind(easy_matches, missing)

# now get the NCBI taxon IDs using 'taxizedb'
ncbi.names <- name2taxid(all_matches$scientificName, db = "ncbi", out_type = "summary")

# aaaaand add in the NCBI names to the dataframe
ncbi_tax <- NULL
for (jj in 1:nrow(all_matches)){
  curr.match <- all_matches[jj,]
  if (curr.match$scientificName %in% ncbi.names$name){ncbi_tax <- append(ncbi_tax, ncbi.names[which(ncbi.names$name == curr.match$scientificName),"id"])}
  else {ncbi_tax <- append(ncbi_tax, NA)}
}
ncbi_tax <- as.vector(unlist(ncbi_tax))
all_matches$taxonID <- ncbi_tax

# make our data look like the BPA metadata sheet
metadata <- data.frame(extraction_id = all_matches$Extraction, # this is a new column and should be removed before submitting to BPA
                       sample_id = "",
                       specimen_id = paste(all_matches$institutionCode, all_matches$catalogNumber),
                       specimen_id_description = all_matches$collectionName,
                       tissue_number = unlist(all_matches$catalogNumber),
                       voucher_or_tissue = unlist(all_matches$catalogNumber),
                       institution_name = all_matches$institutionName,
                       tissue_collection = all_matches$collectionName,
                       sample_custodian = "",
                       access_rights = "No restrictions",
                       tissue_type = "",
                       tissue_preservation = "ethanol",
                       sample_quality = "unknown",
                       taxon_id = all_matches$taxonID,
                       phylum = all_matches$phylum,
                       class = all_matches$class, 
                       order = all_matches$order, 
                       family = all_matches$family,
                       genus = all_matches$genus,
                       species = sapply(all_matches$species, function(x) str_split(x, " ")[[1]][2]),
                       lineage = all_matches$scientificName, # this is a new column and should be removed before submitting to BPA
                       subspecies = "",
                       common_name = all_matches$common_name_and_lsid,
                       identified_by = all_matches$identifiedBy,
                       collection_date = all_matches$eventDate,
                       collector = all_matches$recordedBy,
                       collection_method = "",
                       wild_captive = "wild",
                       source_population = NA,
                       country = all_matches$country,
                       state_or_region = all_matches$stateProvince,
                       location_text = all_matches$locality,
                       habitat = all_matches$habitat,
                       decimal_latitude = all_matches$decimalLatitude,
                       decimal_longitude = all_matches$decimalLongitude,
                       coord_uncertainty = all_matches$coordinateUncertaintyInMeters,
                       genotypic_sex = "not determined",
                       phenotypic_sex = all_matches$sex,
                       method_of_determination = "",
                       certainty = "",
                       life_stage = all_matches$lifeStage,
                       birth_date = NA,
                       death_date = NA,
                       associated_media = "",
                       ancillary_notes = NA,
                       barcode_id = "",
                       ALA_specimen = all_matches$recordID,
                       prior_genetics = "",
                       taxonomic_group = "Squamata",
                       type_status = all_matches$typeStatus,
                       Material_extraction_type = "DNA",
                       Material_extraction_date = "",
                       Material_extracted_by = "",
                       Material_extraction_method = "Salt extraction",
                       Material_conc_ng_ul = "",
                       Genome_sample = "",
                       Genome_status = "",
                       Phylogenomic_sample = "",
                       Phylogenomic_status = "",
                       Conservation_sample = "",
                       Conservation_status = "")

# write your metadata to a csv file, and add it to the sheet here: 
# https://docs.google.com/spreadsheets/d/1MG3DbOEgEtww2S1Y6XBmxg_0YWS1mh1Faf1jsl4QLi8/edit?ts=60cae464#gid=0
write.csv(metadata, file="~/Desktop/Pygopodidae_Carphodactylidae_Metadata.csv", row.names = F)
