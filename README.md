# DeST_OT
Developmental Spatiotemporal Optimal Transport

![figure_1](https://github.com/raphael-group/DeST_OT/blob/main/fig1.png)

A method for alignment of spatially resolved transcriptomics time-series. 


There are three main functions:
1. `src/destot/DESTOT/align`: Given a pair of ST slices from two developmental timepoints, find a spatiotemporal alignment matrix Pi and a growth vextor xi.
2. `src/destot/metrics/growth_distortion_metric`: Given a pair of ST slices, their spatiotemporal alignment matrix, and the inferred growth vector, calculcate the growth distortion metric as in Eq. 9 of the paper.
3. `src.destot/metrics/migration_metric`: Given a pair of ST slices and their spatiotemporal alignment matric, calculate the migration metric as in Eq. 11 of the paper.

## Installation
We will soon make DeST-OT available on PyPi. In the mean time, you can download the repository and call the functions directly.

## Contact
If you encounter any problem running the software, please contact Xinhao Liu at xl5434@princeton.edu or Peter Halmos at ph3641@princeton.edu
