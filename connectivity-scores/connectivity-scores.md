Connectivity scores
================
Pieter Moris
21 februari 2016

Unraveling the connection between chemical structure and gene expression.
-------------------------------------------------------------------------

The aim of this analysis is to investigate the association between the structure of chemical compounds (or their predicted biological activity) and the effects they induce in a gene expression assay. We will group compounds with a similar structure (*chemical fingerprint*) or predicted bioactivity (*target prediction*) using clustering methods. Subsequently the *gene expression profiles* of the compounds within a cluster will be contrasted with each other using a connectivity score approach based on the work of [Lamb et al. (2006)](#lamb) and [Zhang & Gant (2008)](#zang) on the [Connectivity Map](https://www.broadinstitute.org/cmap/).

Their methods enable us to contrast a number of query expression profiles with a set of reference profiles and assess the *connectivity* between compounds in terms of which genes are regulated up or down. The connection between the query and reference data can be both positive and negative. In the former case, highly up-regulated genes in one set will also be up-regulated in the other set, and vice versa for down-regulated genes. In the latter case, genes with a high expression in one assay, might be strongly down-regulated in the other one. Both situations imply that the two compounds interfere with the same biological processes. If two sets are weakly connected, there is no correlation or overlap between the top genes of either set.

*To do: description and origins of dataset.*

Clustering chemical compounds
-----------------------------

We will employ hierarchical clustering to group the chemical compounds based on their fingerprints or target predictions. In either case we are dealing with a 0/1 matrix where the rows correspond to different compounds. For the chemical fingerprints the columns indicate the presence or absence of specific chemical structures. Together this sequence of ones and zeros describes the entire molecular scaffolding of the compound. In essence we are dealing with a bit string for each compound.

*To do: add description of bioactivity matrix + describe dimensions of data.*

Let's take a look at the chemical fingerprints first. Here is a sneak peak at the data:

    ##                  -2147375257 -2147119955 -2146474760 -2145840573
    ## metformin                  0           0           0           0
    ## phenformin                 0           0           0           0
    ## phenyl biguanide           0           0           0           0
    ## estradiol                  0           0           0           0
    ## dexamethasone              0           0           0           0
    ## verapamil                  0           0           0           0

To properly cluster objects based on binary attributes, we need to define an adequate measure of similarity The *Tanimoto coefficient* (sometimes called the *Jaccard coefficient*) is often used for this purpose in cheminformatics ([MacCuish and MacCuish, 2011](#chem)). It contrasts two objects in terms of the number of common attributes, *c*, with the number of attributes unique to either object, *a* and *b*. \[i^n X_i\]

Promising reading materials
---------------------------

-   [Use of chemical similarity in drug discovery](http://mcc.irb.hr/mcc_04/presentations/butina_d_mcc04_2.pdf)
-   [ChemmineR](https://www.bioconductor.org/packages/3.3/bioc/vignettes/ChemmineR/inst/doc/ChemmineR.html#introduction) and [tutorial](http://chemmine.ucr.edu/help/):
    -   contains modules for similarity searching of chemical compounds (Tanimoto coefficient)
    -   clustering of compounds by structural and physicochemical similarities is a powerful approach for correlating structural features of compounds with their activities
-   [textbook on clustering in drug discovery](https://books.google.be/books?id=ZDDNBQAAQBAJ&pg=PA41&lpg=PA41&dq=tanimoto+clustering&source=bl&ots=vsLen2ZmS5&sig=N16boAKkB5NWLJjeteC4shM6Brc&hl=en&sa=X&ved=0ahUKEwiaiZvsx4vLAhWBfxoKHf-9DrwQ6AEISzAH#v=onepage&q=tanimoto%20clustering&f=false)
-   Some cross-validated links:
    -   [frequent item set better than clustering?](https://stats.stackexchange.com/questions/86318/clustering-a-binary-matrix)
    -   [ordinal vs nominal binary](https://stats.stackexchange.com/questions/116856/hierarchical-or-two-step-cluster-analysis-for-binary-data?rq=1)
    -   [hierarchical clustering](https://stats.stackexchange.com/questions/2717/clustering-with-a-distance-matrix?rq=1)
    -   [more clustering options for boolean](https://stats.stackexchange.com/questions/70113/cluster-large-boolean-dataset?rq=1)

References
----------

<a name="lamb"></a>- Lamb, J. 2006. "The Connectivity Map: Using Gene-Expression Signatures To Connect Small Molecules, Genes, And Disease". Science 313 (5795): 1929-1935. <doi:10.1126/science.1132939>.

<a name="chem"></a>- MacCuish, John D, and Norah E MacCuish. 2011. Clustering In Bioinformatics And Drug Discovery. Boca Raton: Taylor & Francis.

<a name="zang"></a> - Zhang, Shu-Dong, and Timothy W Gant. 2008. "A Simple And Robust Method For Connecting Small-Molecule Drugs Using Gene-Expression Signatures". BMC Bioinformatics 9 (1): 258. <doi:10.1186/1471-2105-9-258>.
