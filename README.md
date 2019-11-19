# README

## The ABBA-BABA test

The original ABBA-BABA test for gene flow, conceived by Nick Patterson (Green et al. 2011), is based based on a hypothetical phylogeny with one distant outgroup (H4) and three ingroup samples (H1,H2,H3,H4).

          /\
         /  \
        /\   \
       /  \   \
      /\   \   \
     /  \   \   \
    H1  H2  H3  H4

    A   B   B   A
    B   A   B   A

This approach uses haplotype sequencing data and counts the occurence of either of two genotype patterns "A-B-B-A" and "B-A-B-A" in the four samples, which are not conistent with the tree topology (see below). In the case of incomplete lineage sorting, both ABBA and BABA should occur at equal proportions.

<img src=/images/DG.png height="55" />

Thus count(ABBA)-count(BABA), scaled by the total number of counts (*D*<sub>G</sub>; see above) should equal 0. However, if there is geneflow from H3 to H2, we expect an excess of ABBA, resulting in negative values of *D*<sub>G</sub>.

Equivalent to Green *et al.*’s (2010) test, Durand *et al.* (2011; *D*<sub>D</sub>) and Soraggi *et al.* (2018; *D*<sub>S</sub>) recently proposed using allele frequencies instead of counts. Here, *p*<sub>ij</sub> is the frequency of the derived allele B at the *i*<sup>th</sup> locus in the *j*<sup>th</sup> population in the generic phylogeny (P1, P2, P3 or P4).

<img src=/images/DD.png height="50" />

<img src=/images/DS.png height="50" />

## ABBA-BABA-4-AF

The python script "ABBA-BABA-4-AF.py" is a new implementation of the ABBA-BABA test for gene flow based on the latter two approaches and calculates both genome-wide and window-wise estimates of *D* based on a given (fixed) number of SNPs per window. A 2D matrix of allele frequencies is used as the input where columns represent samples and rows SNP positions. Following Green *et al.* (2011), the software also calculates *z*-scores based on jackknifing using a blocked, even *m*-delete jackknife procedure, where *m* is the group size of observations removed from the sample for jackknifing (Busing et al. 1999). See our corresponding manuscript "*Genomic signals of admixture as well as reinforcement between two sister species of European sepsid flies in sympatry vs. parapatry*" (Giessen *et al.* 2020) for more details on the approach and differences to other methods

## Workflow

The ususal starting point for the analysis of geneflow with ABBA-BABA-4-AF is the generation of an Allele frequency matrix file. This can be done, for example, by converting a SNP file in the SYNC file format (Kofler *et al.* 2011) to major allele frequencies using the provided python script "

## References

@busingDeletemJackknifeUnequal1999

@greenDraftSequenceNeandertal2010

@durandTestingAncientAdmixture2011

@koflerPoPoolation2IdentifyingDifferentiation2011

@soraggiPowerfulInferenceDStatistic2018
