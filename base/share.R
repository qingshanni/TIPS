
##### public functions 
readGMT <- function (file) 
{
  if (!grepl("\\.gmt$", file)[1]) {
    stop("Pathway information must be a .gmt file")
  }
  geneSetDB = readLines(file)
  geneSetDB = strsplit(geneSetDB, "\t")
  names(geneSetDB) = sapply(geneSetDB, "[", 1)
  geneSetDB = lapply(geneSetDB, "[", -1:-2)
  geneSetDB = lapply(geneSetDB, function(x) {
    x[which(x != "")]
  })
  return(geneSetDB)
}



getInputValue <- function(v,v.default){
  if(is.null(v)){
    flag <- v.default
  }else{
    flag <- v
  }
}

sigmoid <- function(pst, params,n=100) {
  
  min_max <- range(pst)
  p <- min_max[1] + (min_max[2] - min_max[1])*(0:n)/n
  mu0 <- params[1] 
  k <- params[2] 
  t_0 <- params[3]
  mu <- 2 * mu0 / (1 + exp(-k*(p - t_0)))
  return(data.frame(x=p,y=mu))
}

getBreaks <- function(v,n){
  min_max <- range(v)
  n <- n-1
  p <- min_max[1] + (min_max[2] - min_max[1])*(0:n)/n
  p
}


