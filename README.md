# ABBA-BABA-4-AF

The python script "ABBA-BABA-4-AF" is a new implementation of the ABBA-BABA test for gene flow based on the approaches of Durand et al. (2011) and Soraggi et al. (2018).

The original ABBA-BABA test for gene flow, conceived by Nick Patterson (Green et al. 2011), is based based on a hypothetical phylogeny with one distant outgroup (H4) and three ingroup samples (H1,H2,H3,H4). This approach uses haplotype sequencing data and counts the occurence of either of two genotype patterns "A-B-B-A" and "B-A-B-A" in the four samples, which are not conistent with the tree topology (see below). In the case of incomplete lineage sorting, both ABBA and BABA should occur at equal proportions. Thus count(ABBA)-count(BABA), scaled by the total number of counts (ABBA+BABA) should equal 0. 

![Pattersons D](/images/DG.png)
However,  under gene-flow from H3 to H2, the





          /\
         /  \
        /\   \
       /  \   \
      /\   \   \
     /  \   \   \
    H1  H2  H3  H4
    
    A   B   B   A
    B   A   B   A
