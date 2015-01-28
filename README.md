# Morphological diversity in tenrecs (Afrosoricida, Tenrecidae): Comparing tenrec skull diversity to their closest relatives
[Sive Finlay](http://sivefinlay.com/) and [Natalie Cooper](https://www.tcd.ie/Zoology/research/ncooper/nataliecooper.php).

This repository contains all of the data and code for the manuscript

## Manuscript
The latest version of the manuscript is available in the Manuscript folder.

To compile the paper:

```
make -C Manuscript
```

## Data
The data is available in the data folder (sub-divided into separate data sets for skulls in dorsal (skdors), ventral (skvent) and lateral (sklat) views.
There are three files for each analysis: a TPS file with the landmarks and curves, an NTS sliders file identifying the curves and an excel file with taxonomic information for each image.
Photographs of some of the skulls are available on [Figshare](http://figshare.com/authors/Sive_Finlay/414410). Museum copyright restrictions prevented public sharing of other images but they are available on request.

## Analysis
There are two R scripts in the analysis folder.

diversity_twofamily contains all of the code to clean the data, conduct geometric morphometrics analyses and estimate morphological diversity in each Family. There are (commented) options to choose which analysis to run (all species/ with a subsample of the Microgale tenrecs) and with which data set.

diversity_PCAplots_3panel contains the code used to create the PCA results figure in the paper
