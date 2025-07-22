# Codes for realization of sub-ångström resolution ptychography in a scanning electron microscope at 20 keV

Arthur M. Blackburn, University of Victoria, 2025

## 1. Introduction

The codes given here produce the key results of work related to the realization of sub-ångström resolution ptychography in a scanning electron microscope at 20 keV, which has recently been submitted for publication[1]. The codes provides a wrapper for and modifies a version of Ptychoshelves_EM  [2, 3], which in turn is a modified subset of a prior release of PtychoShelves[4]. 

To allow some traceability back to the Ptychoshelves_EM  [2, 3] code, the first commit to the repository here is a copy of the code hosted on Zenodo [3], with some superfluous backup (\*.bak, \*.asv) files removed, and some – but not all - Matlab based Unix / Windows mixed formatting resolved. Modifications we have made to the Ptychoshelves_EM code can be seen by viewing the repository changelog. 

The wrapper and modifications provided herewith and called ‘PtychoRunner’, facilitate using PtychoShelves_EM and PtychoShelves in a batch run environment, where all the parameters controlling the reconstructions are contained within a \*.xls spreadsheet file. These files can be produced by applications such as MS Excel, LibreOffice, OpenOffice Calc, or similar. PtychoRunner sets the parameters for reconstructions in a single controlling spreadsheet and allows the output of one ptychographic reconstruction run (or recipe) to be the initial object and probe model to another subsequent recipe. The aim of the wrapper is to facilitate tracking and control of reconstructions and datasets, in small to medium-sized, manual user-led investigations of reconstruction parameter spaces and experiments. It is particularly helpful for running PtychoShelves on batch scheduled systems, such as Slurm, such as used on Digital Research Alliance of Canada computing facilities. Alternatively, automated parameter space searches can be performed using entries from the \*.xls sheet as a baseline and using the wrapper to selectively vary a subset of the parameters.

PtychoRunner has been developed mainly for internal research group usage at UVic, and in a similar style to PtychoShelves is documented in an ad-hoc basis within the code. As with PtychoShelves, ability to trace parameters in the provided code, is helpful or required to understand the function of the control parameter. In due course the documentation provided here may be improved or this code may be forked or branched for future research work. Searching the Matlab \*.m files using text based search tools, such grep, Agent Ransack or similar, to find the parameter names (see the provided example \*.xls sheet column headers) in the code can assist in the process of understanding the function of the parameters. Many of these parameters are briefly explained in the `conv_ptychoset_cSAXS_MSlc` method of `ptychoset_mult.m`.

In addition to the above, the code in PtychoRunner also:

- Calculates the pincushion distortion present in diffraction data by fitting either electron powder diffraction data or the Fourier transform of the ptychographic reconstruction to model electron diffraction data.
- Applies pincushion distortion correction on diffraction data.
- Determines the Fourier Ring Correlation (FRC) between images of samples that have an incomplete or sparse distribution of structured particles. The method involves sub-sampling the larger field of view, then applying Otsu’s [5] method to determine a threshold that minimize the inter-class variance on the mean FRC scores for the sub-region. Regions dominated by the amorphous support membrane, and/or possibly mobile or time varying regions that typically have an FRC score below this threshold, are then excluded from an average FRC characteristic.
- Computes the Spectral Signal to Noise Ratio (SSNR) and Phase Contrast Transfer Function (PCTF) of a reconstruction if the ground truth is known, such as is the case the non-experimental computationally determined diffraction data.

To illustrate the above, and reproduce the key results of the manuscript, example scripts are given within this repository that: 

1. Demonstrate creating a simple reconstruction from synthetically created data. This also provides a quick test to check the installation and functioning of the core code.
2. Creates the ptychographic reconstructions given in the manuscript, specifically using:  
   1. Gold (Au) particles on thin amorphous carbon (aC) experimental diffraction data. 
   2. Two subsets of the Au particles on aC experimental data. These sets have a reduced sampling density (over the image plane), that the above source dataset, but have larger lateral span in the image plane.
   3. Gold island on few layer MoS2 experimental data.
   4. Various synthetic datasets used to create Extended Data Figure 8 in the related manuscript.
   5. Simulated diffraction data produced from a hexagonal frustrum gold island on 5 layers of MoS2, where the entire system is tilted slightly to best match our experimentally determined tilt angle.
3. Perform a fit of the reconstructed Au on aC data to model polycrystalline Au diffraction data and perform a subsequent distortion-correction of the diffraction data. Note: Interested readers may also want to investigate our more recent code and related publications for correcting diffracted data, that does not require model fitting [https://github.com/mrfitzpa/emicroml].
4. Determine the Fourier Ring Correlation (FRC) characteristic between the Au on aC reconstructions, using the method described in manuscript.
5. Determine the Spectral Signal to Noise Ratio (SSNR) and Phase Contrast Transfer Function (PCTF) of the reconstructed simulated datasets, as also given in Extended Data Figure 8 in the related manuscript.

Related to the code provide here, the source data required to produce the reconstructions etc, are provided at a separate repository [https://doi.org/10.20383/103.01357, available after review embargo period]. This dataset also contains the output reconstructions and some intermediate results.

## 2. Description of Example Scripts and Data Repository

Each of the example scripts outlined in the introduction is now described in turn, after first describing the general setup procedure for the software. 

### 2.1    General Setup of Paths

The first step is to setup your local paths. The main directory that your matlab installation must be able to find is the 'Runner' directory of the repository. This can be set through various means:

#### 2.1.1 Matlab Startup File

Place a startup.m file in Matlab's 'userpath' folder. The userpath folder is top of Matlab's search path (aside from the Matlab's current working directory, which is always the very first in the search order), so any files placed in the userpath - such as `startup.m` - override other same-named files placed elsewhere on the search path. To find which folder is the userpath use the command `userpath` at the Matlab command prompt. In the created `startup.m` file you can then add a line that points to the contents of the 'Runner' directory of the repository. For example, add the following line to the `startup.m` file:

    addpath('/projects/my_repos/PtychoRunner/Runner');

#### 2.1.2 Path Definition

Another way to modify the path above is to add the \\Runner path to the Matlab pathdef file. Often a copy of pathdef is also stored in the userpath directory,

    home_directory/Documents/MATLAB/pathdef.m 

which overrides the Matlab’s default pathdef which on a Windows platform is typically installed at 

    C:\Program Files\MATLAB\R2xxx\toolbox\local\pathdef.m.

Use the Matlab command to `which pathdef` to see which pathdef is in use.

However, usually this file is written using the Matlab GUI. Mathworks do not recommend manually editing this file using a text editor, but opening it is a quick way to see what paths are being used. Additionally, permissions on the installation directory of Matlab may prevent you from directly editing a version of pathdef.m that lies in your system related directories, so the version in your userpath (if it exists) will be the editable one that matters.

#### 2.1.3 Edit Runner/SetupPathsN as necessary for your platform

With the \Runner path added to Matlab and active, check and edit if necessary the contents of the `SetupPathsN.m` file which is within the \Runner directory. 

Running the Matlab command

    exec_params = SetupPathsN(true);

should further alter Matlab's path so that the leading entries point to various points within repository. SetupPathsN is a convenience function to set up run parameters appropriately for the platform you are using.

Make sure that the `exec_param` structure entries ending in `_path` point to locations within the repository for example:

> `ptycho_matlab_path: 'repos\PtychoRunner\PtychoShelvesMS\ptycho'`
> `cSAXS_matlab_path: ‘repos\PtychoRunner\PtychoShelvesMS'`

It is also worthwhile setting up a useful default location for where dataset information will be located. For example, this might be:

    datasets_path: 'projects\Ptycho\ExampleDataSets\records'

### 2.2 Description of Example Scripts

The \TestScripts directory contains example scripts pertaining to the numbered items in the introduction. In addition to the numbered items described in the introduction, a ‘zeroth’ example script is also given, which describes our initial pre-processing of the data.

#### 2.2.1 Data Preparation and Preprocessing (Part_0_Prepare_Data.m)

This script demonstrates the pre-processing of the 'raw' data obtained in HDF5 format from Hitachi High-Tech's Azorus software into a Matlab v7.3 \*.mat file (also HDF5 based) format that is used in the subsequent demonstration files. Briefly, the script:

- Opens the raw data
- Recentres it, using either manual selection of the direct beam, or using some ‘automatic’ options such as centre of mass.
- Forms subsets of the collected diffraction data. Two types of subsets are formed: one where the number of diffraction sampling points is reduced, but the sampling density on the image plane remains the same; and another where the number of diffraction sampling points and the sampling density are both halved.
- Shows how distortion correction is applied to the MoS2 diffraction data set. This uses distortion coefficients determined from `Part_3_DistortionFit_and_Correct.m` (see 2.1.3), which determines the distortion present in the Au particles in amorphous carbon data.

#### 2.1.2 Create of Reconstructions (Part_1and2_Reconstructions.m)

This script demonstrates how to run reconstructions based on parameters and datasets described in an excel file. The script also allows the creation of the all the reconstructions referred to and analyzed within the text.

The script contains code blocks that call the `RunReconF` function. This function initiates a reconstruction using a reconstruction recipe and dataset that is specified within the excel file given by ‘ExcelFileName’ function argument. More details on thus use of RunReconF are found within that function. Briefly though, the recipe and data to be reconstructed, is specified by the `exec_char` function argument. This is a single character (or combination of up to 4 characters) that should match the `do_recon` and `exec` columns of the Dataset and Recipes sheets of the provided excel sheet.

Blocks are provided to initiate runs for the key reconstructions, specifically the blocks provide reconstructions of:

1. Synthetic data, providing a quick test run to see if things are working okay. This test run should take 1 - 2 minutes, or less, to run and requires < 2 GB of GPU memory.
2. Synthetic data at conditions very similar to those used in Au on MoS2 reconstruction.
3. Synthetic data at conditions very similar to those used in the Au on amorphous carbon (aC) reconstruction.
4. Experimental data from Au on Amorphous carbon (aC) - with no distortion correction.
5. Experimental data from Au on Amorphous carbon (aC) - after first round of distortion correction.
6. Experimental data from Au on Amorphous carbon (aC) - after second (and final) round of distortion correction.
7. Experimental data from Gold Islands on MoS2 (using distortion corrected data)
8. Simulated from a hexagonal frustrum gold island on 5 layers of MoS2, where the entire system is tilted slightly to best match our experimentally determined tilt angle.

Some of these blocks will need to be un-commented in the Matlab file in order to be active and run. Furthermore, the blocks and parameters might need to be adjusted for your architecture. The first 'quick test' block, requires 1 GPU, whereas the next block require a computer with 3 GPUs (or 3 virtual GPUs on the system), in order to run 3 reconstructions in parallel.

To modify the second block to run on just one GPU change `exec_char`, `rst`,  to `exec_char`, `r`. The characters `r,s,t` in this example relate to datasets and recipes in the Excel sheet. If `rst` is specified, datasets `r, s, t` with related recipes `r, s, and t` run in parallel. If just `r` is specified, then only dataset `r` is processed with the related recipes `r`.

The synthetic dataset (#1 above) can be preformed on a GPU having 2 GB or less of dedicated memory. Systems without a GPU can also be used (but have not been thoroughly tested), but obviously run slower. 
The output from a run produces two output files.

1. A summary monochrome tif format image file. This file allows quick inspection of the features of the produced reconstruction, and is divided into 6 regions:
   
   - Left hand side, upper row: The amplitude of the reconstructed probe wave function.
   - Left hand side, lower row: The phase of the reconstructed probe wave. (Note the phase is not unwrapped).
   - Middle, upper row: The phase of reconstructed object. In the case of a multi-slice model the presentation is by default the exit wave phase produced by propagating an initial plane-wave through the slice of the model. As with the probe, the phase is not unwrapped in this presentation.
   - Middle, lower row: The Fourier transform (FT) of the phase given in the row above. If there are significant (unphysical) phase wrapping artefacts in the above image, then this will have an effect on the FT.
   - Right hand side, upper row: The amplitude of the formed reconstruction. In the case of a multi-slice reconstruction the single plane presented is the same as used for the phase image given in the middle upper row.
   - Right hand side, lower row: The FT of the amplitude presented above.
     Note, however, that the amplitude images have 1% of the extremal pixel values saturated by default in order to improve image contrast in the majority of the image; see savesummaryimage3.m, so the unprocessed reconstruction data stored in the \*.mat file should be used in preference to the image for further numerical analysis. 

2. A matlab m v7.3 \*.mat file (which is a compressed HDF5 based format), that contains the full precision data of the completed reconstruction including the full multi-slice object model and potential multi-mode probe model. Numerous other fields are also within the structure, that describe the parameters used by PtychoShelves and its modifications to define the reconstruction. (Note however that there may be some redundancy in fields of the Matlab structure, and some fields might be set to default values that do not actively specify the reconstruction parameters used).This file also contains a record of the progression of the Fourier error metric of the reconstruction with iteration number to provide a trace of the convergence (or potentially divergence) of the reconstruction with iterations.

Images of the reconstructed probe modes and object slices at intermediate stages are also produced as the reconstruction is progressing if the ‘save_results_every’ column of the recipe is set to a number less than the total number of iterations in the recipes step. These can be used to monitor to reconstruction progress or see how the reconstruction developed over the iterations in the reconstruction.

#### 2.1.3. Part_3_DistortionFit_and_Correct.m

This script provides an example of how to determine and correct for pincushion distortion from a reconstruction of randomly oriented gold particles. This was as applied to the data as described in the manuscript, and in the comments provided within this script. 

#### 2.1.4. Part_4_FRC_Characteristics.m

This script re-runs the calculations used to produce the FRC characteristics given in the manuscript.

#### 2.1.5. Part_5_PCTF_SNR_on_SyntheticData.m

This script determines the SNR and PCTF from the ptychographic data reconstructed on synthetic test objects, where the simulated data has the same or very similar reconstructed pixel size and similar (or larger) 20 keV electron beam illumination half-angle as in the experimental data.

## 3. Description of Dataset Repository Contents

The source data and key reconstructions produced from the codes provided are available in a publicly accessible dataset hosted by the Canadian Federated Research Data Repository [https://doi.org/10.20383/103.01357, available after review embargo period]. The description of the data contained in the dataset is provided in the file `TestScripts\Related_Dataset_Readme.txt`.

## References

1. Blackburn, A. M., Cordoba, C., Fitzpatrick, M. R. C. and Mcleod, R. A., 2025,  Sub-ångström resolution ptychography in a scanning electron microscope at 20 keV, Submitted, Under Review.
2. Chen, Z., Jiang, Y., Shao, Y.-T., Holtz, M. E., Odstrčil, M., Guizar-Sicairos, M., Hanke, I., Ganschow, S., Schlom, D. G. and Muller, D. A., 2021,  Electron ptychography achieves atomic-resolution limits set by lattice vibrations, Science 372, 826.
3. Chen, Z. J., Yi; Muller, David A.; Odstrčil, Michal (2021) PtychoShelves_EM, source code for multislice electron ptychography, (Zenodo  https://doi.org/10.5281/zenodo.4659690 ).
4. Wakonig, K., Stadler, H.-C., Odstrcil, M., Tsai, E. H. R., Diaz, A., Holler, M., Usov, I., Raabe, J., Menzel, A. and Guizar-Sicairos, M., 2020,  PtychoShelves, a versatile high-level framework for high-performance analysis of ptychographic dataThis article will form part of a virtual special issue of the journal on ptychography software and technical developments, Journal of Applied Crystallography 53, 574-586.
5. Otsu, N., 1979,  A Threshold Selection Method from Gray-Level Histograms, IEEE Transactions on Systems, Man, and Cybernetics 9, 62-66.
