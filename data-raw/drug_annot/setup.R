library(data.table)
library(dplyr)
library(tidyr)
library(drugseqr)

pkgconfig::set_config("dplyr::na_matches" = "never")

# function definitions ----

# non UTF-8 encoded character cause DT alerts on filtering
remove_non_utf8 <- function(df) {
  df[] <- lapply(df, function(x) {
    Encoding(x) <- "UTF-8"
    iconv(x, "UTF-8", "UTF-8",sub='')
  })

  return(df)
}

summarise_func <- function(x) {

  unqx <- na.omit(x)
  unqx <- unlist(strsplit(unqx, '|', fixed = TRUE))
  unqx <- unique(unqx)

  # keep as NA if they all are
  if (!length(unqx)) return(NA_character_)

  # collapse distinct non-NA entries
  return(paste(unqx, collapse = ' | '))
}

rename_cols <- function(annot) {
  annot <- rename(annot,
                  `Clinical Phase` = clinical_phase,
                  MOA = moa,
                  Target = target,
                  `Disease Area` = disease_area,
                  Indication = indication,
                  Vendor = vendor,
                  `Catalog #` = catalog_no,
                  `Vendor Name` = vendor_name)

  return(annot)
}

destructure_title <- function(pdata, drop = TRUE, ...) {
  title_split <- strsplit(pdata$title, '_')

  pdata <- tibble::add_column(pdata,
                              'Compound'  = sapply(title_split, `[`, 1),
                              'Cell Line' = sapply(title_split, `[`, 2),
                              'Dose'      = sapply(title_split, `[`, 3),
                              'Duration'  = sapply(title_split, `[`, 4), ...)

  if (drop) pdata$title <- NULL
  return(pdata)
}

setup_annot <- function(annot, study) {

  increasing_phases <- c('Withdrawn', 'Preclinical', 'Phase 1', 'Phase 2', 'Phase 2 | GRAS', 'Phase 3', 'Launched', 'Launched | GRAS')

  pdata <- readRDS(paste0('inst/extdata/', study, '_pdata.rds'))
  pdata <- destructure_title(pdata, drop=FALSE, .after=1)

  # pubchem_cid and pert_iname together have the most matches
  sum(pdata$Compound %in% tolower(annot$pert_iname))
  sum(pdata$`Pubchem CID` %in% annot$pubchem_cid)
  sum(pdata$`Pubchem CID` %in% annot$pubchem_cid | pdata$Compound %in% tolower(annot$pert_iname))

  # first get matches on cid
  cat('Matching on cids ... \n')
  annot_cids <-  pdata %>%
    select(title, Compound, 'Pubchem CID') %>%
    inner_join(annot, by = c('Pubchem CID' = 'pubchem_cid'))

  # where there wasn't a cid match try on pert_iname
  cat('Matching on Broad internal names ... \n')
  study_annot <- pdata %>%
    select(title, 'Pubchem CID', 'Compound') %>%
    anti_join(annot, by = c('Pubchem CID' = 'pubchem_cid')) %>%
    left_join(annot, by = c('Compound' = 'pert_iname')) %>%
    bind_rows(annot_cids)


  # fill in na values with non-na value by cid and name
  cat('Filling in NA by CIDs and compound names ... \n')
  study_annot <- study_annot %>%
    group_by(`Pubchem CID`) %>%
    fill(everything()) %>%
    fill(everything(), .direction = 'up') %>%
    ungroup() %>%
    group_by(Compound) %>%
    fill(everything()) %>%
    fill(everything(), .direction = 'up') %>%
    ungroup()

  # use most advanced phase if same title has multiple phases
  cat('Selecting most advanced clinical phase for each title ... \n')
  phase <- study_annot %>%
    select(title, clinical_phase) %>%
    mutate(clinical_phase = factor(clinical_phase, increasing_phases, ordered = TRUE)) %>%
    group_by(title) %>%
    summarise(clinical_phase = max(clinical_phase))

  # merge rows that have the same title
  cat('Merging rows with the same title ... \n')
  study_annot <- study_annot %>%
    select(-clinical_phase) %>%
    group_by(title) %>%
    summarise_all(summarise_func) %>%
    left_join(phase, by = 'title') %>%
    select(title, `Pubchem CID`, clinical_phase, everything()) %>%
    arrange(match(title, pdata$title))

  # check that study_annot and pdata are in same order
  stopifnot(all.equal(study_annot$title, pdata$title))

  # save as seperate annotation, removing what pdata already has
  study_annot <- study_annot %>% select(-c(title, `Pubchem CID`, pert_iname, pubchem_cid, Compound)) %>% rename_cols()
  saveRDS(study_annot, paste0('inst/extdata/', study, '_annot.rds'))

  return(study_annot)
}

# setup BRH data ----

drugs <- fread('data-raw/drug_annot/repurposing_drugs_20180907.txt', skip = 'pert_iname')
samples <- fread('data-raw/drug_annot/repurposing_samples_20180907.txt', skip = 'pert_iname')

# fix up some imperfections
table(samples$pert_iname %in% drugs$pert_iname)
# FALSE  TRUE
#     2  10145

samples$pert_iname[!samples$pert_iname %in% drugs$pert_iname]
# "golgicide-A" "YM-298198-desmethyl"

grep('golgicide', drugs$pert_iname, value = TRUE)
grep('golgicide', samples$pert_iname, value = TRUE)
samples[samples$pert_iname == 'golgicide-A', 'pert_iname'] <- 'golgicide-a'

# merge based on pert_iname
brh_annot <- left_join(drugs, samples, by = 'pert_iname')
brh_annot[brh_annot == ''] <- NA

# keep what is relevant for now
brh_annot <- select(brh_annot, -c(broad_id, qc_incompatible, purity, expected_mass, smiles, InChIKey, deprecated_broad_id))
brh_annot <- unique(brh_annot)
brh_annot$pubchem_cid <- as.character(brh_annot$pubchem_cid)

# use most advanced phase
brh_annot$clinical_phase <- gsub('^Phase \\d/', '', brh_annot$clinical_phase)
unique(brh_annot$clinical_phase)

brh_annot <- remove_non_utf8(brh_annot)


# add SIDER, DrugBank, Wikipedia, and GRAS ----

sider <- readRDS('data-raw/drug_annot/pug_view/sider.rds')
sider <- tibble(pubchem_cid = names(sider), SIDER = sider) %>%
  mutate(SIDER = ifelse(sider, pubchem_cid, NA))

pug_annot <- readRDS('data-raw/drug_annot/pug_view/pug_annot.rds')
pug_annot <- rename(pug_annot, DrugBank = drugbank, Wikipedia = wikipedia)

# get matches on pubchem cid
annot <- brh_annot %>%
  as_tibble() %>%
  select(-pert_iname) %>%
  left_join(sider, by = 'pubchem_cid') %>%
  inner_join(pug_annot, by = 'pubchem_cid') %>%
  mutate(clinical_phase = ifelse(gras, paste0(clinical_phase, ' | GRAS'), clinical_phase)) %>%
  select(-gras) %>%
  select(clinical_phase, DrugBank, SIDER, Wikipedia, everything())

# get those that didn't match on pubchem cid and match by pert_iname
annot <- brh_annot %>%
  as_tibble() %>%
  anti_join(pug_annot, by = 'pubchem_cid') %>%
  select(-pubchem_cid) %>%
  full_join(pug_annot, by = c('pert_iname')) %>%  # need full_join for subsequent fill
  select(-gras) %>%
  bind_rows(annot)

# fill in non-NA by cid then by pert_iname
annot <- annot %>%
  group_by(pubchem_cid) %>%
  fill(everything()) %>%
  fill(everything(), .direction = 'up') %>%
  ungroup() %>%
  group_by(pert_iname) %>%
  fill(everything()) %>%
  fill(everything(), .direction = 'up') %>%
  ungroup() %>%
  distinct()

# CMAP02/L1000 setup ----

cmap_annot <- setup_annot(annot, 'CMAP02')
l1000_annot <- setup_annot(annot, 'L1000')

