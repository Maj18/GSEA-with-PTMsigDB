# How to run GSEA on phosphosite data using PTMsigDB



## Download [PTMsigDB](https://proteomics.broadapps.org/ptmsigdb/](https://github.com/broadinstitute/ssGSEA2.0/tree/master/db/ptmsigdb/v2.0.0/all)



## Prepare the database

```
ptmsigdb = GSEABase::getGmt("../data/db_ptm.sig.db.all.v2.0.0/ptm.sig.db.all.flanking.human.v2.0.0.gmt")
ptmsigdb_KINASE0 = ptmsigdb[grepl("KINASE", names(ptmsigdb))]
ptmsigdb_KINASE0 = lapply(ptmsigdb_KINASE0, function(K) {
  data.frame(Kinase=K@setName, Phosphosite=K@geneIds)
}) %>% Reduce(rbind, .)
ptmsigdb_KINASE2 = as.data.frame(stringr::str_split_fixed(ptmsigdb_KINASE0$Phosphosite,";",2)) %>%
  setNames(c("Phosphosite", "Regulation"))
dim(ptmsigdb_KINASE2)
dim(ptmsigdb_KINASE0)
ptmsigdb_KINASE = data.frame(Kinase=ptmsigdb_KINASE0$Kinase,
                             Phosphosite=ptmsigdb_KINASE2$Phosphosite,
                             Regulation=ptmsigdb_KINASE2$Regulation)
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
                                         TERM2GENE = ptmsigdb_KINASE)
```
