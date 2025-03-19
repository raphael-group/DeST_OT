# DeST_OT
**De**velopmental **S**patio**T**emporal **O**ptimal **T**ransport

![figure_1](https://github.com/raphael-group/DeST_OT/blob/main/fig1.png)

A method for aligning spatially resolved transcriptomics time-series. 


There are four main functions:
1. `src/destot/DESTOT/align`: Given a pair of ST slices from two developmental timepoints, infer a spatiotemporal alignment matrix Pi and a growth vector xi. As discussed in the paper, there are two settings we recommend for this function: the default setting for growth-rates ($\alpha = 0.2$, $\beta = 0.5$, $\gamma = 50$, $\epsilon = 0.1$), and the setting ($\alpha = 0.99$, $\beta = 0.6$, $\gamma = 1$, $\epsilon = 0.1$) for spatiotemporal alignments with differing geometries (e.g. with capture-frame effects). Larger entropy $\epsilon$ is more accurate for growth-rates, and smaller (down to 1e-4) is more appropriate for one-to-one alignments.
2. `src/destot/DESTOT/xi_to_growth_rate`: Given a growth vector xi, convert the values in the growth vector to a per-spot growth rate $J$ given the start and end timepoints.
3. `src/destot/metrics/growth_distortion_metric`: Given a pair of ST slices, their spatiotemporal alignment matrix Pi, and the inferred growth vector xi, calculcate the growth distortion metric as in Eq. 9 of the paper.
4. `src.destot/metrics/migration_metric`: Given a pair of ST slices and their spatiotemporal alignment matrix Pi, calculate the migration metric as in Eq. 11 of the paper.

## Installation
We will soon make DeST-OT available on PyPi. In the mean time, you can download the repository and call the functions directly.

## Contact
If you encounter any problem running the software, please contact Xinhao Liu at xl5434@princeton.edu or Peter Halmos at ph3641@princeton.edu

## Reference
Halmos, P., Liu, X., Gold, J., Chen, F., Ding, L., and Raphael, B. J. DeST-OT: Alignment of Spatiotemporal Transcriptomics Data. Cell Systems, January 2025. ISSN 2405-4712. doi: 10.1016/j.cels.2024.12.001. URL http://dx.doi.org/10.1016/j.cels.2024.12.001.

The paper is available here: <[DeST-OT: Alignment of Spatiotemporal Transcriptomics Data.](https://www.cell.com/cell-systems/fulltext/S2405-4712(24)00365-X)>, and a Zenodo registered DOI is available in the link below

[![DOI](https://zenodo.org/badge/714103322.svg)](https://zenodo.org/doi/10.5281/zenodo.13769696)

