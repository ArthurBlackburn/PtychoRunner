PtychoShelves_EM

This source code is developed based on an older version of PtychoShelves package from Paul Scherrer Institut (PSI). It was modified for multislice electron ptychographic reconstruction at Cornell University, USA. The related paper (or subsequent journal version) should be cited whenever this code, dataset or their derivatives are used: 
Chen, Z., et al., Electron ptychography achieves atomic-resolution limits set by lattice vibrations. arXiv: 2101.00465.
Wakonig K. et al., PtychoShelves, a versatile high-level framework for high-performance analysis of ptychographic data, J. Appl. Cryst. 53(2) (2020).
Odstrcil M, et al., Iterative least-squares solver for generalized maximum-likelihood ptychography, Optics express. 2018 Feb 5;26(3):3108-23. 

--- 2021/3/27, Dr. Zhen Chen, Dr. Yi Jiang, Prof. David A. Muller

#######################

- The code is streamlined from the original PtychoShelves code for demonstration purpose only. 

- The complete package with many more algorithms and documentation can be found at https://www.psi.ch/en/sls/csaxs/software

- Copyright and license issues should follow the agreements in PSI's codes and/or refer to PtychoShelves website.  

#######################

Short notes for startup:

1. One experimental example dataset is put in ../PtychoShelves_EM/ptycho/exampleData/  (large diffraction data needs to be loaded from PARADIM website: https://data.paradim.org/, check later after our paper is pulished).

2. For an initial test of the reconstruction, change the Matlab current Folder (work path) to ../PtychoShelves_EM/ptycho/, and run the main drive scriptmo.m. 

3. For your own data, prepare a data file 'data_dp.mat' containing diffractions, 'data_position.hdf5' containing probe positions, and initial probe 'probe_initial.mat'. Please see the example script: ../ptycho/utils_EM/prepare_data_electron.m. 

4. look into the script ptychography_demo.m and modify parameters accordingly. More options can be refered to the full documentation for PtychoShelves. 

5. Only Matlab 2018a with the GPU engine has been tested. 

6. Further improvements may be acheived by using more probe modes, more layers, better sampling, or a better initial start of probe and probe positions. But keep in mind that using large dimensions is very compuational demanding in both memory and time.

#########################

Contact information:

Zhen Chen (zhen.chen@conell.edu) or David A. Muller (david.a.muller@cornell.edu)
