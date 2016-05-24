#' This function can be used to visualice the number of GO (Gene Ontology) signals
#' gained and lost between the WGCNA standard network and the new network obtained
#' from the application of k-means search heuristic. It is based on the gProfileR
#' package
#' @param net.label A label to use within the plots, to refer to the GCN
#' @param net The network obtained by the k-means process (WGCNA object)
#' @param partitions The partition file
#' @param plot.file If not NULL, then it should be the pdf file name in which the
#' GO gain/lost will be presented.
#' @param filter The filter used for gProfileR (see details at this package)
#' @param exclude.iea If true, do not use Inferred Electronic Annotations in the
#' enrichment analysis (gProfileR)
#' @param gprof.method The default is the one implemented by gProfiler (see package
#' for details)
#' @return The results from the analysis for further study
#' @importFrom gProfileR gprofiler

plotGOenrichment <- function(net.label,
                             net,
                             partitions,
                             plot.file=NULL,
                             filter="GO",
                             exclude.iea=T,
                             gprof.method="gSCS"){


  gene.names <- names(net$moduleColors)
  background = gene.names
  modules <- unique(net$moduleColors)
  module.colors = modules

  enrichment <- matrix(nrow=length(partitions),ncol=length(modules))
  new.terms <- matrix(nrow=length(partitions),ncol=length(modules))
  lost.terms <- matrix(nrow=length(partitions),ncol=length(modules))
  highly.signif.terms <- matrix(nrow=length(partitions),ncol=length(modules))

  colnames(enrichment) <- modules
  colnames(new.terms) <- modules
  colnames(lost.terms) <- modules

  enrichment[,] <- 0

  partition.index <- 1

  goresults <- list()

  #Adjust 1st partition to the rest
  names(partitions[[1]]) <- names(partitions[[2]])

  names(partitions[[length(partitions)]]) = names(partitions[[2]])
  partitions.to.plot <- partitions[c(1,length(partitions))]

  for(partition in partitions.to.plot){
    print(paste0("Going for enrichment in iteration ", partition.index))

    all.genes <- NULL
    for(module in modules){

      all.genes[[module]] = names(partition)[partition == which(module.colors == module)]
    }
    go <- gprofiler(query=all.genes,correction_method=gprof.method,exclude_iea=exclude.iea,
                    custom_bg=background,src_filter=filter)

    goresults[[partition.index]] <- go
    #Transforming p values in log10 scale
    en.modules <- unique(go$query.number)
    for(module in en.modules){
      enrichment[partition.index,module] <- sum(-log10(go$p.value[go$query.number == module]))
      if(partition.index == 1){
        new.terms[partition.index,module] <- 0
      }else{
        old.go <- goresults[[partition.index - 1]]
        intersect.go <- intersect(old.go$term.id[old.go$query.number == module],
                                  go$term.id[go$query.number == module])
        new.terms[partition.index,module] <- length(setdiff(go$term.id[go$query.number == module],
                                                            intersect.go))

        lost.terms[partition.index,module] <- length(setdiff(old.go$term.id[old.go$query.number == module],
                                                             intersect.go))
        new.terms[is.na(new.terms[,])] <- 0
        lost.terms[is.na(lost.terms[,])] <- 0
      }
    }
    partition.index <- partition.index + 1
  }

  if(!is.null(plot.file)){
    pdf(plot.file,width=20,height=10)
    old.par <- par()
    par(mfrow=c(1,3))
  }
  vals.to.plot <- enrichment[2,]-enrichment[1,]
  vals.to.plot <-  c(vals.to.plot,sum(enrichment[2,])-sum(enrichment[1,]))
  order.to.plot = order(vals.to.plot,decreasing=TRUE)
  mod.cols <- c(modules,"gold")
  barplot(vals.to.plot[order.to.plot],col=mod.cols[order.to.plot],
          main=paste0(net.label, ": GO gain after post-processing WGNCA"),xlab="Network modules",
          ylab="sum(-log10(pval))" )

  vals.to.plot <- new.terms[2,]-new.terms[1,]
  vals.to.plot <-  c(vals.to.plot,sum(new.terms[2,])-sum(new.terms[1,]))
  order.to.plot = order(vals.to.plot,decreasing=TRUE)
  barplot(vals.to.plot[order.to.plot],col=mod.cols[order.to.plot],
          main=paste0(net.label,": new terms by module"),xlab="Network modules",ylab="# terms")

  vals.to.plot <- lost.terms[2,]-lost.terms[1,]
  vals.to.plot <-  c(vals.to.plot,sum(lost.terms[2,])-sum(lost.terms[1,]))
  order.to.plot = order(vals.to.plot,decreasing=TRUE)

  barplot(vals.to.plot[order.to.plot],col=mod.cols[order.to.plot],
          main=paste0(net.label,": lost terms by module"),xlab="Network modules",ylab="# terms")

  if(!is.null(plot.file)){
    par(old.par)
    dev.off()
  }

  return(list(enrichment=enrichment,goresults=goresults))
}
