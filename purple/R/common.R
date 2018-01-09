left_join<-function(left, right, by="sampleId") {
  tmp = merge(x = left, y = right, by=by, all.x=TRUE)
  return (tmp)
}

slookup<-function(df, lookup, on="sampleId") {
  return (sapply(df[[on]], function(x) {lookup[match(x, lookup[[on]]), -c(1)] }))
}

#this is a slow function - don't use!
dflookup<-function(df, lookup, on="sampleId") {
  return(as.data.frame(t(sapply(df[[on]], function(x) {lookup[match(x, lookup[[on]]), -c(1)] }))))
}
