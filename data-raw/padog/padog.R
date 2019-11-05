gslist <- crossmeta::gslist
gs.names <- crossmeta::gs.names

saveRDS(gslist, 'data-raw/padog/gslist.rds')
saveRDS(gs.names, 'data-raw/padog/gs.names.rds')


# get for GO terms
ont <- select(GO.db, columns = c('GOID', 'ONTOLOGY'), keys = keys(GO.db))
ont <- ont[ont$ONTOLOGY == 'BP', ]
gslist.go <- as.list(org.Hs.egGO2ALLEGS)
gslist.go <- gslist.go[names(gslist.go) %in% ont$GOID]
gs.names.go <- AnnotationDbi::Term(names(gslist.go))

hs <- readRDS("data-raw/homologene/hs.rds")

# subset to hs enids and use same symbols as cmap_es and l1000_es
gslist.go <- lapply(gslist.go, function(x) {x[x %in% hs$ENTREZID]})
gslist.go <- lapply(gslist.go, function(x) {names(x) <- hs[x]$SYMBOL_9606; x})

gslist.go <- gslist.go[sapply(gslist.go, length) >= 3]
gs.names.go <- gs.names.go[names(gslist.go)]


saveRDS(gslist.go, 'data-raw/padog/gslist.go.rds')
saveRDS(gs.names.go, 'data-raw/padog/gs.names.go.rds')
