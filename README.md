# Reconstruction of consensus metabolic models and community-dependent gap-filling using COMMIT for communites sampled from _Arabidopsis thaliana_.

## Description
This repository contains all files used for iterative gap-filling and subsequent analysis
of 432 consensus genome-scale metabolic models for OTUs sampled from _A. thaliana_ roots, leaves and surrounding soil.
The consensus models were generated by merging draft models from four different approaches:
- [KBase](https://www.kbase.us/) (Arkin et al. 2018)
- [RAVEN 2.0](https://github.com/SysBioChalmers/RAVEN) (Wang et al. 2018)
- [CarveMe](https://github.com/cdanielmachado/carveme) (Machado et al. 2018)
- [AuReMe](http://aureme.genouest.org/) (Aite et al. 2018) / [Pathway Tools](http://pathwaytools.com/) (Karp et al. 2016)


## Requirements for COMMIT
- Matlab (tested with version R2017b)
- [CPLEX solver](https://www.ibm.com/analytics/cplex-optimizer) (tested with v12.9)
- R packages: pheatmap, wesanderson, scales, plotrix, igraph, ape
- [COBRA toolbox v3.0](https://github.com/opencobra/cobratoolbox)
- [HMMER](http://hmmer.org/download.html) (tested with v3.2.1)

## Reference