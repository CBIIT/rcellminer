searchIdsOutput <- DT::renderDataTable({
  drugAnnot <- as(featureData(getAct(rcellminerData::drugData)), "data.frame")[, c("NSC", "NAME", "MOA")]

  dtColnames <- c("NSC", "Name", "MOA")
  tmp <- DT::datatable(drugAnnot, rownames=FALSE, colnames=dtColnames, filter='none', style='bootstrap', selection='none', options=list(pageLength=10))
})
