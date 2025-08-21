#!/usr/bin/env Rscript
# r_sequence_collation.R
# Usage:
# Rscript r_sequence_collation.R --input parsed/cox1_parsed.csv --family_lookup data/family_lookup.csv --out results/cox1_by_marker.xlsx

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(stringr)
})

option_list = list(
  make_option(c("-i","--input"), type="character", help="Input parsed CSV/XLSX"),
  make_option(c("-f","--family_lookup"), type="character", default=NULL, help="CSV mapping Genus->Family (optional)"),
  make_option(c("-o","--out"), type="character", default="results/collated.xlsx")
)
opt = parse_args(OptionParser(option_list=option_list))

input_file = opt$input
family_lookup_file = opt$family_lookup
out_file = opt$out

if(grepl("\\.xlsx?$", input_file)){
  df = read_excel(input_file)
} else {
  df = read.csv(input_file, stringsAsFactors = FALSE)
}

# Ensure sequence col exists
if(!"RNA Sequence" %in% colnames(df)){
  stop("Input must contain 'RNA Sequence' column")
}

# Helper functions
extract_main <- function(org){
  if(is.na(org)) return(NA)
  org = str_trim(as.character(org))
  m = str_match(org, "^\\s*([A-Za-z][A-Za-z\\-]+)")
  if(!is.na(m[1,2])) return(tolower(m[1,2]))
  return(tolower(org))
}
extract_species <- function(org){
  if(is.na(org)) return(NA)
  parts = str_split(org, "\\s+")[[1]]
  if(length(parts) >= 2) return(paste(parts[-1], collapse = " "))
  return(NA)
}

df$Organism = as.character(df$Organism)
df$Main_Organism = sapply(df$Organism, extract_main)
df$Species = sapply(df$Organism, extract_species)
df$`Length (nt)` = nchar(gsub("-", "", df$`RNA Sequence`))

# Family lookup
family_map = data.frame()
if(!is.null(family_lookup_file) && file.exists(family_lookup_file)){
  family_map = read.csv(family_lookup_file, stringsAsFactors = FALSE)
  # expect columns: Genus, Family
  colnames(family_map) = tolower(colnames(family_map))
  if("genus" %in% names(family_map) && "family" %in% names(family_map)){
    lookup = setNames(tolower(family_map$family), tolower(family_map$genus))
    df$Family = lookup[df$Main_Organism]
  }
}
df$Family[is.na(df$Family)] = "unknown"

# cp and f-h mapping if exists
ecol = data.frame()
if(file.exists("data/family_ecology_lookup.csv")){
  ecol = read.csv("data/family_ecology_lookup.csv", stringsAsFactors = FALSE)
  ecol$family = tolower(ecol$Family)
  ecol_map_cp = setNames(ecol$cp, ecol$family)
  ecol_map_fh = setNames(ecol$`f-h`, ecol$family)
  df$cp = ecol_map_cp[df$Family]
  df$`f-h` = ecol_map_fh[df$Family]
}
df$cp[is.na(df$cp)] = "Unknown"
df$`f-h`[is.na(df$`f-h`)] = "Unknown"

# Save annotated
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
write.xlsx(df, out_file)

# produce simple summary
summary_df = df %>%
  group_by(Family) %>%
  summarise(NumSeq = n(), NumGenera = n_distinct(Main_Organism), NumSpecies = n_distinct(Species)) %>%
  arrange(desc(NumSeq))
write.xlsx(summary_df, file = sub("\\.xlsx$", "_summary.xlsx", out_file))

cat("Collation complete. Output:", out_file, "\n")
