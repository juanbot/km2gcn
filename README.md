# km2gcn
Optimization process of WGCNA hierarchical clustering with k-means

This package is an additional step for refining the gene clusters obtained by WGCNA from a TOM (Topological Overlap Measure). By default, WGCNA uses hierarchical clustering, using complete linkage and a distance matrix based on the TOM. Once WGCNA finishes, it generated a network with module eigengenes and a partition, represented as a vector of genes, in the names() attribute, and the modue colors of each cluster as the vector content.

This package starts with such object, and applying a k-means heuristic, improves the clusters in many directions:
-it increases the eigengene as a proxy by getting more genes in each cluster whose MM is the highest in that module.
-it increases the GO enrichment of the modules (it uses gProfileR) to generate a GO functional description of module function
-it increases module preservation.

A paper which describes the approach is on the way
