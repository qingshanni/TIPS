cat('\f')
rm(list=ls())

#dd <- as.matrix(read.csv( "D:\\work\\Zihan\\资料\\error05082020\\tcellraw.csv",row.names = 1))




load('tt.RData')

aa <- sdata@meta.data


gene_num_cutoff <- 0  ###
gene_set.new    <- lapply(gset, function(genes){intersect(genes, rownames(sde)) })
gene_set.new.n  <- lapply(gene_set.new,length)
gene_set.new    <- gene_set.new[gene_set.new.n >= gene_num_cutoff]

pathway_sde <- lapply(gene_set.new,function(genes){ 
  s <- sde[genes,] 
  rownames(s) <- s$gene
  s
})
names(pathway_sde) <- names(gene_set.new)







gene_set <-  ptime[['gene_set']]   ## use gene set filter in Pseudotime
pathway_sde <- lapply(gene_set,function(genes){ 
  s <- sde[genes,] 
  rownames(s) <- s$gene
  s
  })
names(pathway_sde) <- names(gene_set)
sdes <- pathway_sde

sig  <- lapply(sdes, function(sde) sde[(sde$qval <0.05), ])  #removes non-significant genes
N    <- sapply(sig, function(s) dim(s)[1])    #gene number
t0   <- sapply(sig, function(s) median(s$t0))    #median timepoint for switch
k    <- sapply(sig, function(s) median(abs(s$k))) #median switch intensity



mx     <- gdata_expr1
#mx     <- mx[rowSums(mx)>0,]
ptime  <- pData(phenoData(scData1))$Pseudotime
gene_set <-  PseudotimeListAll1$gene_set   ## use gene set filter in Pseudotime
#sde      <- lapply(gene_set, function(genes) switchde(mx[genes, ], ptime)) #returns a list of dfs

sde <- list()
idx <- c()
n <- length(gene_set)
k <- 0
for(i in 1:n){
  tryCatch({
    tmp <- switchde(mx[gene_set[[i]], ], ptime)
    k <- k+1
    sde[[k]] <- tmp
    idx <- c(idx,i)
    },error=function(e){})
}

names(sde) <- names(gene_set)[idx]

sde1  <- sde[[1]]
genes <- sde1$gene[1:3]
mx.new<- mx[genes,]
dm    <- dim(mx.new)
v     <- as.vector(t(mx.new))

dd    <- data.frame(Expression = v, Pseudotime = rep(ptime,dm[1]),gene = rep(genes,each= dm[2] ))
min_max <- range(ptime)
data.line <- NULL

for(gene in genes){
  df <- sigmoid(ptime,params = extract_pars(sde1, gene))
  data.line <- rbind(data.line,data.frame(df,gene = rep(gene, dim(df)[1])))
}

ggplot() + theme_bw() +
  geom_point(data=dd, aes(x = Pseudotime, y = Expression,color = gene,fill=gene), shape = 21) +  
  geom_line(data=data.line, aes(x = x, y = y,color = gene),size=1)






