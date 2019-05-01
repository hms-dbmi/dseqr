library(data.table)
library(dplyr)

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


# CMAP02 pdata
cmap_pdata_orig <- readRDS('inst/extdata/CMAP02_pdata.rds')
cmap_pdata_orig <- destructure_title(cmap_pdata_orig, drop=FALSE, .after=1)

# pubchem CID seems to have the most matches
table(unique(cmap_pdata_orig$Compound) %in% annot$pert_iname)
table(unique(cmap_pdata_orig$`Pubchem CID`) %in% annot$pubchem_cid)

# use cmap_pdata compound names
cmap_pdata <- left_join(cmap_pdata_orig, annot, by = c('Pubchem CID' = 'pubchem_cid'))

# merge rows that have the same treatment
cmap_pdata <- cmap_pdata %>%
  group_by(title) %>%
  summarise_all(function(x)  {
    unqx <- unique(na.omit(x))

    # keep as NA is they all are
    if (!length(unqx)) return(NA_character_)

    # collapse distinct non-NA entries
    return(paste(unqx, collapse = ' // '))
  }) %>%
  arrange(match(title, cmap_pdata_orig$title))


all.equal(cmap_pdata$title, cmap_pdata_orig$title)
