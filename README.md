# Steering-Floquet
This repository contains Jupyter notebooks for reproducing the results in our article "Steering spin fluctuations in lattice systems via two-tone Floquet engineering" published in Physica Scripta https://iopscience.iop.org/article/10.1088/1402-4896/ad9d85. Our work considers a one-dimensional spin-1/2 lattice with periodically modulated bonds using parametric resonances. Quantum dynamics is studied using exact diagonalization (ED) up to L=10 spins with the Quspin package https://quspin.github.io/QuSpin/. We also provide a Python script for computing the stroboscopic dynamics.  

Summary of results:

1. We report on the control of spin pair fluctuations using two-tone Floquet engineering.
2. The stroboscopic dynamics generated from distributed spin exchange modulations lead to spin pair fluctuations reaching quasi-maximally correlated states and a subharmonic response in local observables, breaking the discrete-time translational symmetry.
3. We present a protocol to control the interacting many-body dynamics, producing spatial and temporal localization of correlated spin pairs via dynamically breaking correlated spin pairs from the edges towards the center of the lattice.
4. Our result reveals how spin fluctuations distribute in a heterogeneous lattice depending on parametric resonances. 

# Requirements

The library is entirely developed in Python 3 using Numpy and Scipy, and a standard Jupyter environment. We recommend using Anaconda or Miniconda, although any other distribution of Python should suffice. Also, we recommend using the tqdm package for progress bar visualization. 
