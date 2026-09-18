## =============================================================================
## 00 -- Prepare combined phenotype data (Aberystwyth UK + Braunschweig/JKI Germany)
## -----------------------------------------------------------------------------
## Selects the 125 genotypes used in the original sparse-testing analysis from
## the full ABR33 diversity panel (954 genotypes), combines them with the
## Braunschweig/JKI phenotypes, and produces one long-format phenotype table
## spanning all 5 environments (Aberystwyth 2015-2016, Braunschweig 2014-2016),
## with population group (DAPC-derived) attached to every record.
##
## Inputs (raw, unmodified):
##   ../../Original phenotypes.xlsx   (sheet "ABR33Data_ForNelson") -- UK, 954 genotypes
##   ../../JKI_Ger.csv                                              -- Germany
##   ../data/popgroups_125.csv -- the 125-genotype list and its DAPC PopGroup
##       assignment, derived from marker data by 00a_derive_population_groups.R
##       (run that script first if this file does not yet exist)
##
## Output:
##   data/combined_phenotypes_long.csv
## =============================================================================

rm(list = ls())
suppressMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

setwd("/Users/nel6/Desktop/Documents/Miscathus_others/Sparse_testing/revision/Miscanthus_SparseTesting_Pipeline/scripts")

RAW_UK_XLSX   <- "../../Original phenotypes.xlsx"
RAW_JKI_CSV   <- "../../JKI_Ger.csv"
GENOLIST_CSV  <- "../data/popgroups_125.csv"
OUT_CSV       <- "../data/combined_phenotypes_long.csv"

## -----------------------------------------------------------------------------
## 1. Authoritative 125-genotype list + population group (marker-derived,
##    see 00a_derive_population_groups.R)
## -----------------------------------------------------------------------------
if (!file.exists(GENOLIST_CSV)) {
  stop("Run 00a_derive_population_groups.R first -- ", GENOLIST_CSV, " does not exist.")
}
genotypes_125 <- read.csv(GENOLIST_CSV, stringsAsFactors = FALSE) |>
  distinct(Genotype, PopGroup) |>
  arrange(Genotype)
stopifnot(nrow(genotypes_125) == 125)
cat("Loaded 125-genotype list (marker-derived). PopGroup counts:\n")
print(table(genotypes_125$PopGroup))

## -----------------------------------------------------------------------------
## 2. UK / Aberystwyth phenotypes (2015, 2016)
## -----------------------------------------------------------------------------
uk_raw <- read_excel(RAW_UK_XLSX, sheet = "ABR33Data_ForNelson")
uk_raw <- uk_raw |> filter(!is.na(geno))

## The raw file contains 13 exact-duplicate UID rows (same plant record
## entered twice) across the full 954-genotype panel; drop them before
## selecting our 125 genotypes, or the affected plots are double-counted.
n_before <- nrow(uk_raw)
uk_raw <- uk_raw |> distinct(UID, .keep_all = TRUE)
cat("UK: removed", n_before - nrow(uk_raw), "exact-duplicate UID rows\n")

uk_sel <- uk_raw |>
  filter(geno %in% genotypes_125$Genotype) |>
  select(Genotype = geno, block, row, col,
         dry_matter_plant_15, fresh_weight_plant_15, moisture_content_15,
         dry_matter_plant_16, fresh_weight_plant_16, moisture_content_16)

cat("\nUK: matched", n_distinct(uk_sel$Genotype), "of 125 genotypes\n")

uk_long <- uk_sel |>
  pivot_longer(
    cols = c(dry_matter_plant_15, fresh_weight_plant_15, moisture_content_15,
             dry_matter_plant_16, fresh_weight_plant_16, moisture_content_16),
    names_to = c("trait", "yr"),
    names_pattern = "(.*)_(\\d+)$",
    values_to = "value"
  ) |>
  mutate(
    Year = paste0("20", yr),
    Environment = paste0("Aber-", Year),
    Rep = NA_character_,      # UK design uses block/row/col; no separate replicate term
    Block = as.character(block),
    Row = as.character(row),
    Col = as.character(col),
    value = suppressWarnings(as.numeric(value))
  ) |>
  select(Genotype, Environment, Year, Rep, Block, Row, Col, trait, value)

## -----------------------------------------------------------------------------
## 3. Braunschweig / JKI phenotypes (2014, 2015, 2016)
## -----------------------------------------------------------------------------
## Genotype naming differs between sites: JKI uses "Mb 1037#52", UK/125-list
## uses "Mb-1037_52". Normalize JKI names to the UK convention before matching.
jki_raw <- read.csv(RAW_JKI_CSV, stringsAsFactors = FALSE)
jki_raw$Genotype_norm <- jki_raw$Genotype |>
  str_replace_all(" ", "-") |>
  str_replace_all("#", "_")

## Exact-duplicate records (same UID, date, and trait values entered twice)
n_before_jki <- nrow(jki_raw)
jki_raw <- jki_raw |> distinct(UID, Pheno_date, Year, Dry_matter_plant, Fresh_weight_plant, Moisture_content, .keep_all = TRUE)
cat("JKI: removed", n_before_jki - nrow(jki_raw), "exact-duplicate rows\n")

jki_sel <- jki_raw |>
  filter(Genotype_norm %in% genotypes_125$Genotype) |>
  filter(!is.na(Dry_matter_plant) | !is.na(Fresh_weight_plant) | !is.na(Moisture_content))

cat("\nJKI: matched", n_distinct(jki_sel$Genotype_norm), "of 125 genotypes\n")
cat("JKI records by year (after filtering to 125 genotypes):\n")
print(table(jki_sel$Year))

## Parse "JKI 15 AM-R2 C1" -> Row = 2, Col = 1. Row numbering reveals a
## genuine two-block replicated design (not a single-plant design as
## previously assumed): the same genotype at the same Col reappears at
## Row+18 in ~90% of cases across all three years (verified: 2014
## 192/209, 2015 185/198, 2016 95/107 matching genotype x Col pairs),
## i.e. Row 1-18 = Block 1, Row 19-35 = Block 2.
jki_sel <- jki_sel |>
  mutate(
    Row = str_extract(Location, "R\\d+") |> str_remove("R") |> as.integer(),
    Col = str_extract(Location, "C\\d+") |> str_remove("C"),
    Block = ifelse(Row <= 18, "1", "2")
  )

jki_long <- jki_sel |>
  select(Genotype = Genotype_norm, Year, Block, Row, Col,
         dry_matter_plant = Dry_matter_plant,
         fresh_weight_plant = Fresh_weight_plant,
         moisture_content = Moisture_content) |>
  mutate(Year = as.character(Year), Row = as.character(Row)) |>
  pivot_longer(cols = c(dry_matter_plant, fresh_weight_plant, moisture_content),
               names_to = "trait", values_to = "value") |>
  mutate(
    Environment = paste0("Brau-", Year),
    Rep = NA_character_
  ) |>
  select(Genotype, Environment, Year, Rep, Block, Row, Col, trait, value)

## -----------------------------------------------------------------------------
## 4. Combine, attach population group, save
## -----------------------------------------------------------------------------
combined <- bind_rows(uk_long, jki_long) |>
  left_join(genotypes_125, by = "Genotype") |>
  filter(!is.na(value))

## Consistent trait abbreviations used throughout the pipeline and all
## downstream outputs: DM = dry matter, FW = fresh weight, MC = moisture
## content.
trait_abbrev <- c(dry_matter_plant = "DM", fresh_weight_plant = "FW", moisture_content = "MC")
stopifnot(all(combined$trait %in% names(trait_abbrev)))
combined$trait <- trait_abbrev[combined$trait]

cat("\n=== Combined dataset summary ===\n")
cat("n genotypes:", n_distinct(combined$Genotype), "\n")
cat("Environments:\n"); print(table(combined$Environment))
cat("Records per trait x environment:\n")
print(table(combined$trait, combined$Environment))
cat("\nAny genotype missing PopGroup after join?\n")
print(sum(is.na(combined$PopGroup)))

write.csv(combined, OUT_CSV, row.names = FALSE)
cat("\nSaved:", OUT_CSV, "\n")
