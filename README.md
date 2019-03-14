# km2gcn
Optimization process of WGCNA hierarchical clustering with k-means

Before you keep on reading, please consider moving from this package to CoExpNets. If it a km2gcn with steroids. Works better, has more abundant documentation and you can incorporate precomputed co-expression networks in your analyses.
Check <http://github.com/juanbot/CoExpNets>.


km2gcn is not maintained anymore.


This package is an additional step for refining the gene clusters obtained by WGCNA from a TOM (Topological Overlap Measure). By default, WGCNA uses hierarchical clustering, using complete linkage and a distance matrix based on the TOM. Once WGCNA finishes, it generated a network with module eigengenes and a partition, represented as a vector of genes, in the names() attribute, and the modue colors of each cluster as the vector content.

This package starts with such object, and applying a k-means heuristic, improves the clusters in many directions:

-it increases the eigengene as a proxy by getting more genes in each cluster whose MM is the highest in that module.

-it increases the GO enrichment of the modules (it uses gProfileR) to generate a GO functional description of module function

-it increases module preservation.

A paper which describes the approach is on the way

To install from R console, issue these commands

*****

library(devtools)

install_github(repo="juanbot/km2gcn/km2gcn")

*****

Alternatively you can download the source tarball from

https://github.com/juanbot/km2gcn/blob/master/km2gcn_0.1.0.tar.gz

and install from the source. 

And here is an example

*****

library(km2gcn)

data(km2gcndata)

net = applykM2WGCNA(net.label="dummy", net.file=km2gcndata$net, expr.data=km2gcndata$expr, job.path="~/tmp/", meg=0) 

*****
  
For more information on how the algorithm works, check the paper here. And if you want to reference as, use this paper as the reference.

https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-017-0420-6

*****
