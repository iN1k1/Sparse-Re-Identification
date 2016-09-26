This package provides you an updated version of the MATLAB code for following paper: 
Martinel, N., & Micheloni, C. (2014). Sparse Matching of Random Patches for Person Re-Identification. In Proceedings of International Conference on Distributed Smart Cameras, 1–6. doi: 10.1145/2659021.2659034

Please notice that this is not the exact replica of the code we have used to obtain the results resported in the paper.
It has been adapted to make it readable and easy to use with different datasets and features.
As a result, re-identification performance may differ from the ones reported in the paper.

## USAGE:

### Simple Demo
```MATLAB
 startup;
```
Initialize all the directiories which are needed to run the main algorithm

```MATLAB
results = main();
```
Reproduces a single experiment using the WARD dataset. Results, like the CMC and the normalized Area Under Curve are stored in the results structure.

### Settings
If you want to play with the approach parameters, please refer to the 
```MATLAB
init_parameters.m
```
file.

## COMPILE:
Please note that some libraries contain mex-files that needs to be compiled for your machine. We provide a limited set of binary within the package. We have not yet provided a script to compile all the dependences.

## CITATION:
If you use the code contained in this package we appreciate if you'll cite our work. 
BIBTEX:
@inproceedings{Martinel2014a,
author = {Martinel, Niki and Micheloni, Christian},
booktitle = {International Conference on Distributed Smart Cameras},
doi = {10.1145/2659021.2659034},
isbn = {9781450329255},
pages = {1--6},
publisher = {ACM Press},
title = {{Sparse Matching of Random Patches for Person Re-Identification}},
url = {http://dl.acm.org/citation.cfm?doid=2659021.2659034},
year = {2014}
}
