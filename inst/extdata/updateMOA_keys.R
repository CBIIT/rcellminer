# tmpEnv <- new.env()
# 
# load("data/Drug_MOA_Key.RData", envir = tmpEnv)
# dim(tmpEnv$Drug_MOA_Key) # 149   4

Drug_MOA_Key = read.delim("./inst/extdata/moa230_description_cdb.txt",stringsAsFactors = F,row.names=1)
dim(Drug_MOA_Key) # 230   4

# write.table(newMOAkey,"./inst/extdata/Drug_MOA_Key.txt",sep="\t",row.names=F)
# stopifnot(identical(colnames(newMOAkey),colnames(tmpEnv$Drug_MOA_Key)))
# 
save(Drug_MOA_Key, file = "data/Drug_MOA_Key.RData")
