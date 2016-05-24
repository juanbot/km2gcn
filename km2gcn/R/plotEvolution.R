#' This function can be used to visualice the results of a k-means heuristic applied
#' to a WGCNA network. It can be called internally by \code{applykM2WGCNA} or
#' independently.
#' @param net.label A label to use within the plots, to refer to the GCN
#' @param The network obtained by the k-means process (WGCNA object)
#' @param module.colors Final partition (the default will be usually fine)
#' @param tom The TOM (Topology Overlap Matrix) created by WGCNA, or reproduced
#' from the combination of expression data + beta
#' @param beta The soft thresholding used for the GCN in WGCNA
#' @param expr.data The data frame with gene expression (genes in columns, samples in rows)
#' @param net.type  The net type as used in WGCNA
#' @param plot.evolution.file The pdf generated with all the plots
#' @param modules Just if we wanted to focus on specific modules, a list of strings
#' @param go.file One important evaluation of the optimization is the number of Gene Ontology
#' significant terms we get from the modules. If go.file is not NULL, a gProfiler based
#' analysis result will be saved on this file as a rds object. Results will be plotted within
#' the general plot.
#' @return None
#' @export

plotEvolution <- function(net.label,net,
                          module.colors=unique(net$moduleColors),
                          partitions,
                          tom=NULL,
                          beta,
                          expr.data,
                          net.type,
                          plot.evolution.file,
                          modules=NULL,
                          go.file){

  if(is.null(tom)){
    cat("We have to create a TOM matrix from the WGCNA network, with beta",beta,
        "and net type",net.type,"\n")
    #Now we create the TOM
    adjacency = adjacency(expr.data, power = beta, type = net.type )
    tom.matrix = TOMsimilarity(adjacency)
    cat("The tom matrix has been generated\n")
  }else if(typeof(tom) == "character"){
    cat("TOM matrix loaded from file",tom,"\n")
    tom.matrix = readRDS(tom)
  }else{
    cat("TOM matrix passed as argument\n")
    tom.matrix = tom
  }

  pdf(plot.evolution.file,width=15,height=12)
  par(mfrow=c(4,3))
  g.changes = NULL
  n = length(partitions)
  for(index in 1:(n-1)){
    g.changes = c(g.changes,sum(partitions[[index]] != partitions[[index + 1]]))
  }
  plot(1:(length(partitions)-1),g.changes,
       main=paste0(net.label," exchanged genes"),xlab="iteration")

  if(is.null(modules)){
    modules <- 1:length(module.colors)
  }

  if(!is.null(go.file))
    saveRDS(plotGOenrichment(net.label,net,partitions,module.colors),go.file)


  #Lets plot module size evolution
  sizes <- matrix(nrow=length(modules),ncol=n)
  sizes <- lapply(partitions,function(x){
    s <- lapply(1:length(module.colors),function(y,x){
      size <- table(x == y)["TRUE"]
    },x=x)
  })
  min.size <- min(unlist(lapply(partitions,function(x){ return(min(table(x))) })))
  max.size <- max(unlist(lapply(partitions,function(x){ return(max(table(x))) })))

  plot(1:n,lapply(sizes,function(x){ x[[1]]}),
       main=paste0(net.label, " module sizes (",length(sizes)," modules)"),
       col=module.colors[modules[1]],ylim=c(min.size,max.size),t="l")
  for(index in 2:length(modules))
    lines(1:n,lapply(sizes,function(x){ x[[index]]}),col=module.colors[modules[index]])

  if(!is.null(tom.matrix)){
    #Calc adjacency
    adjacency <- matrix(nrow=n,ncol=length(modules))
    max.y <- 0
    min.y <- 1000
    for(module.index in 1:length(modules)){
      for(element in 1:n){
        genes.in.module = partitions[[element]] == modules[module.index]

        #grey will generane NAs. Change to F
        genes.in.module[is.na(genes.in.module)] = F

        if(all(!genes.in.module)){
          print("Module has zero genes ")
        }
        #print(genes.in.module)
        adjacency[element,module.index] = sum(tom.matrix[genes.in.module,genes.in.module])/table(genes.in.module)["TRUE"]
      }
    }

    #Plot adjacency
    max.y <- max(adjacency)
    min.y <- min(adjacency)
    if(is.infinite(max.y))
      max.y = 50
    if(is.infinite(min.y))
      min.y = 0
    cat("Plotting adjacency",c(min.y,max.y),"\n")
    plot(1:n,adjacency[,1],main=paste0(net.label, " adjacency"),ylim=c(min.y,max.y),col=module.colors[modules[1]])
    for(module in modules[-1])
      lines(1:n,adjacency[,module],col=module.colors[module])

    #Plot adjacency entropy
    entropies <- list()
    for(element in 1:n){
      ps <- adjacency[element,]/max(adjacency[element,])
      entropies[[element]] <- -sum(ps * log10(ps))
    }
    cat("Plotting entropy",c(min.y,max.y),"\n")
    plot(1:n,entropies,main=paste0(net.label, " adjacency entropy"),t="l")
  }

  #Plot number of genes exchanged in modules by iteration
  genes.exchanged <- matrix(nrow=(n-1),ncol=length(modules))
  for(module.index in 1:length(modules)){
    for(element in 1:(n - 1)){
      genes.in.module.old <- partitions[[element]] == modules[module.index]
      genes.in.module.new <- partitions[[element + 1]] == modules[module.index]
      genes.exchanged[element,module.index] = table(genes.in.module.old != genes.in.module.new)["TRUE"]
      if(is.na(genes.exchanged[element,module.index]))
        genes.exchanged[element,module.index] <- 0
    }
  }
  max.y <- max(genes.exchanged)
  min.y <- 0
  if(is.infinite(max.y))
    max.y = 4000

  plot(1:(n-1),genes.exchanged[,1],main=paste0(net.label, "gene changes"),
       ylim=c(min.y,max.y),col=module.colors[modules[1]],t="l")
  for(module in modules[-1])
    lines(1:(n-1),genes.exchanged[,module],col=module.colors[module])

  dev.off()

}
