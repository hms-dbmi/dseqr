library(data.table)
library(dplyr)
library(purrr)

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
annot <- left_join(drugs, samples, by = 'pert_iname')
annot[annot == ''] <- NA

# keep what is relevant for now
annot <- select(annot, -c(broad_id, qc_incompatible, purity, expected_mass, smiles, InChIKey, deprecated_broad_id))
annot <- unique(annot)
annot$pubchem_cid <- as.character(annot$pubchem_cid)

# use most advanced phase
annot$clinical_phase <- gsub('^Phase \\d/', '', annot$clinical_phase)
unique(annot$clinical_phase)
increasing_phases <- c('Withdrawn', 'Preclinical', 'Phase 1', 'Phase 2', 'Phase 3', 'Launched')

annot <- remove_non_utf8(annot)

# CMAP02 setup ----
cmap_pdata <- readRDS('inst/extdata/CMAP02_pdata.rds')
cmap_pdata <- destructure_title(cmap_pdata, drop=FALSE, .after=1)

# pubchem CID seems to have the most matches
table(unique(cmap_pdata$Compound) %in% annot$pert_iname)
table(unique(cmap_pdata$`Pubchem CID`) %in% annot$pubchem_cid)

# use cmap_pdata compound names
cmap_annot <-  cmap_pdata %>%
  select(title, 'Pubchem CID') %>%
  left_join(annot, by = c('Pubchem CID' = 'pubchem_cid'))

# use most advanced phase if same CID has multiple phases
cmap_phase <- cmap_annot %>%
  select(title, clinical_phase) %>%
  mutate(clinical_phase = factor(clinical_phase, increasing_phases, ordered = TRUE)) %>%
  group_by(title) %>%
  summarise(clinical_phase = max(clinical_phase))

# merge rows that have the same treatment
cmap_annot <- cmap_annot %>%
  select(-clinical_phase) %>%
  group_by(title) %>%
  summarise_all(summarise_func) %>%
  left_join(cmap_phase, by = 'title') %>%
  select(title, `Pubchem CID`, clinical_phase, everything()) %>%
  arrange(match(title, cmap_pdata$title))

# check that annot and pdata are in same order
all.equal(cmap_annot$title, cmap_pdata$title)

# save as seperate annotation, removing what pdata already has
cmap_annot <- cmap_annot %>% select(-c(title, `Pubchem CID`, pert_iname))
saveRDS(cmap_annot, 'inst/extdata/CMAP02_annot.rds')


# L1000 setup -----

l1000_pdata <- readRDS('inst/extdata/L1000_pdata.rds')
l1000_pdata <- destructure_title(l1000_pdata, drop=FALSE, .after=1)

# pubchem_cid has the most matches
sum(unique(l1000_pdata$Compound) %in% tolower(annot$pert_iname))
sum(unique(l1000_pdata$`Pubchem CID`) %in% annot$pubchem_cid)

# use cmap_pdata compound names
l1000_annot <-  l1000_pdata %>%
  select(title, 'Pubchem CID') %>%
  left_join(annot, by = c('Pubchem CID' = 'pubchem_cid'))

# use most advanced phase if same CID has multiple phases
l1000_phase <- l1000_annot %>%
  select(title, clinical_phase) %>%
  mutate(clinical_phase = factor(clinical_phase, increasing_phases, ordered = TRUE)) %>%
  group_by(title) %>%
  summarise(clinical_phase = max(clinical_phase))


# merge rows that have the same treatment
l1000_annot <- l1000_annot %>%
  select(-clinical_phase) %>%
  group_by(title) %>%
  summarise_all(summarise_func) %>%
  left_join(l1000_phase, by = 'title') %>%
  select(title, `Pubchem CID`, clinical_phase, everything()) %>%
  arrange(match(title, l1000_pdata$title))

# check that annot and pdata are in same order
all.equal(l1000_annot$title, l1000_pdata$title)

# save as seperate annotation, removing what pdata already has
l1000_annot <- l1000_annot %>% select(-c(title, `Pubchem CID`, pert_iname))

# remove non-utf8 and save
saveRDS(l1000_annot, 'inst/extdata/L1000_annot.rds')
