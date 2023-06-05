# sdp_pdes

This repository contains the code for my Spring 2023 EE 364B final project. In it, I implement the method-of-moments for control of non-linear PDEs, first proposed by Korda et al. in 2022 [[1]](#1). My implementation considers the special case of nonlinear HJB equations for a standard optimal control problem with linear-quadratic-Gaussian cost, with the aim of comparing the accuracy and efficiency of the SDP-based approach to HJB equations to a more standard deep learning-based approach such as deep BSDE.

My implementation of the convex optimization-based approach utilizes [GloptiPoly 3](https://homepages.laas.fr/henrion/software/gloptipoly3/), a more general MATLAB framework for polynomial optimization and moment problems. Importantly, GloptiPoly can interface with [YALMIP](https://yalmip.github.io/), a MATLAB toolbox for solving large-scale SDPs using a variety of standard solvers. My approach compares the efficiency of interior-point solvers such as Mosek to that of first-order, ADMM-type solver such as SCS.

## Installation

After cloning this repository, install the [GloptiPoly 3](https://homepages.laas.fr/henrion/software/gloptipoly3/) and [YALMIP](https://yalmip.github.io/download/) MATLAB toolboxes. If you wish to test the Deep BSDE implementation that I've used, originally developed by Han et al. [[2]](#2), you also need to clone their [Deep BSDE](https://github.com/frankhan91/DeepBSDE) repository. Finally, the recommended solver for all SDP relaxations is [SCS](https://www.cvxgrp.org/scs/), which you can download from source or via YALMIP. If you wish to test additional solvers such as Mosek, SDPT3, or SDPNAL, you must download them separately from source or via YALMIP. Once all MATLAB toolboxes and packages have been downloaded, make sure to add all the solvers that you intend on using to your MATLAB path.

## Tests

All test scripts can be found in the `tests` directory. This includes tests for the method-of-moments applied to the HJB equation for an optimal control problem with LQG cost in lower dimensions (i.e., one spatial variable and one temporal variable) *and* scalable versions that can carry out the SDP relaxation procedure in higher dimensions (marked with the tag `scalable`). Note that, due to the sheer size of the moment problems generated for higher dimensions, it may become infeasible to solve the HJB equation in dimensions higher than five. Additionally, I included a simple, one-dimensional test of GloptiPoly 3 in `gloptipoly_test.m`. All test scripts are written in MATLAB. All outputs from the test scripts can be found in the directory `outputs`, and plotting code and final plots are located in `plots`.


## References
<a id="1">[1]</a> 
Milan Korda, D. Henrion, and J.B. Lasserre. 
Moments and convex optimization for analysis and control of nonlinear PDEs. 
*Numerical Control: Part A*, vol. 23, 339-366, 2022.

<a id="2">[2]</a> 
Jiequn Han, A. Jentzen, and E. Weinan. Solving high-dimensional partial differential equations using deep learning.
*Proceedings of the National Academy of Sciences*, 115(34):8505–8510, 2018.
