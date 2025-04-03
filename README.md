# How to run GSEA on phosphosite data using PTMsigDB



## Download [PTMsigDB](https://github.com/broadinstitute/ssGSEA2.0/tree/master/db/ptmsigdb/v2.0.0/all)



## Prepare the database

```
ptmsigdb = GSEABase::getGmt("../data/db_ptm.sig.db.all.v2.0.0/ptm.sig.db.all.flanking.human.v2.0.0.gmt")
# Below we will only extract KINASE from the database, of course you can also extract PERT, DISEASE or PATH as well
ptmsigdb_sub0 = ptmsigdb[grepl("KINASE", names(ptmsigdb))]
ptmsigdb_sub0 = lapply(ptmsigdb_sub0, function(K) {
  data.frame(Set=K@setName, Phosphosite=K@geneIds)
}) %>% Reduce(rbind, .)
ptmsigdb_sub2 = as.data.frame(stringr::str_split_fixed(ptmsigdb_sub0$Phosphosite,";",2)) %>%
  setNames(c("Phosphosite", "Regulation"))
# Only keep the upregulated phosphosites
ptmsigdb_sub2 = ptmsigdb_sub2 %>% filter(Regulation=="up")
ptmsigdb_sub = data.frame(Set=ptmsigdb_sub0$Set,
                             Phosphosite=ptmsigdb_sub2$Phosphosite,
                             Regulation=ptmsigdb_sub2$Regulation)
```


## Prepare the differential phosphosites

```
rslt$PTM_FlankingRegion = paste0(rslt$PTM_FlankingRegion, "-p")
gl = setNames(rslt$logFC, rslt$PTM_FlankingRegion) %>% sort(., decreasing=T)
```



## run GSEA

```
clusterProfiler::GSEA(geneList = gl,
                                         eps=0, 
                                         pvalueCutoff = 1,
                                         pAdjustMethod = "fdr",
                                         TERM2GENE = ptmsigdb_sub)
```
