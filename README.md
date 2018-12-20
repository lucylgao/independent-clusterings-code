# independent-clusterings-code
Code to reproduce simulations and data application from Gao, L.L., Bien, J. and Witten, W. (2018+) Are Clusterings of Multiple Data Views Independent?

## Required R libraries:
mclust
matrixStats
[multiviewtest](https://github.com/lucylgao/multiviewtest)
ggplot2
grid
gridExtra
scales
impute

## Code Dependencies: 
Code for some figures (e.g. Figures 2 and 3 in the main text) depend on the results of auxiliary files. Before running code for the figures, first: 
* Run sim-study-sec4-sec5.R with p = 10 and all combinations of (delta, n, seed, K, sig) in {0, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1} \times {25, 50, 100, 250, 500} \times {1, 2, \dots, 12} \times {3, 6} \times {2.4, 4.8, 9.6}
* Run sim-study-sec4-sec5.R with p = 100 and all combinations of (delta, n, seed, K, sig) in {0, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1} \times {25, 50, 100, 250, 500} \times {1, 2, \dots, 12} \times {3, 6} \times {4.8, 9.6, 19.2}
* Run sim-study-sec4-sec5-collate-results.R with all combinations of (p, K) in {10, 100} \times {3, 6}
* Run chooseK-setting1-sim-code.R, chooseK-setting2-sim-code.R, dense-diagonal-gaussian-sim-code.R, dense-gaussian-sim-code.R, and dense-t-sim-code.R with all combinations of (delta, n, seed) in {0, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1} \times {25, 50, 100, 250, 500} \times {1, 2, \dots, 12}

## Data Dependencies:
Figure 4 and Table 1 in the main text uses data from the Pioneer 100 Wellness Project (Price et. al., 2017), which is publicly available on [dbGAP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/molecular.cgi?study_id=phs001363.v1.p1&phv=301385&phd=&pha=&pht=6518&phvf=&phdf=&phaf=&phtf=&dssp=1&consent=&temp=1). After applying for access to the dbGAP data set and subsequent download of the data, Figure4+Table1.R can be directly applied to participant_data.csv in the main folder, and proteomics_DS.txt, clinical_labs_DS.txt, and metabolites_DS.txt in the dbgap_data subfolder. 
