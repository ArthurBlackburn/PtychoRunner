This README.txt file was generated on 2015-07-18 by Arthur M. Blackburn

--------------------
GENERAL INFORMATION
--------------------

1. Title of Dataset: Electron Diffraction Data and Ptychographic Reconstructions Using a 20 keV Electron Beam on Supported Gold Nanoparticles

2. Author Information
	A. Principal Investigator Contact Information
		Name: Arthur M. BLACKBURN
		Institution:  University of Victoria
		Email: ablackbu@uvic.ca
		Alternative Email: arthur.blackburn@gmail.com

3. Date of data collection
	Experimental Data: Approx 2022-10(Oct)-01 -- 2023-03(Mar)-01
	Synthetic Data and Reconstructions: From above to 2025-04(Apr)
		
4. Geographic location of data collection: 

Advanced Microscopy Facility,
University of Victoria, Victoria, B.C. Canada. 

5. Information about funding sources that supported the collection of the data: 

	This work was part funded by the Natural Sciences and Engineering Research Council of Canada (NSERC), from Collaborative Research and Development grant (CRDPJ 543431), partnering with Hitachi High-Tech Canada. Part of the data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul Scherrer Institute, Switzerland. These codes and supporting simulations were performed using the Advanced Research Computing (ARC) facilities of Digital Research Alliance of Canada. Members of Hitachi High-Tech, Naka, Japan, who developed and integrated the projector lens on the experimental SU9000 instrument are thanked and acknowledged. 

---------------------------
SHARING/ACCESS INFORMATION
---------------------------

1. Licenses/restrictions placed on the data: 

	These data are available under a CC BY 4.0 license <https://creativecommons.org/licenses/by/4.0/> 

2. Links to publications that cite or use the data:
 
	Submitted to Nature Communications, with the title
	"Sub-ångström resolution ptychography in a scanning electron microscope at 20 keV"
	Arthur M. Blackburn 1,2*;
	Cristina Cordoba 1,2;
	Matthew R. Fitzpatrick 1,2;
	Robert A. McLeod 3,4

	1 Department of Physics and Astronomy, University of Victoria, Victoria, BC, Canada
	2 Centre for Advanced Materials and Related Technologies, University of Victoria, Victoria, BC, Canada
	3 Hitachi High-Tech Canada Inc., Toronto, ON, Canada
	4 Department of Materials Science & Engineering, University of Toronto, Toronto, ON, Canada
	
	* Corresponding Author

3. Links/relationships to ancillary data sets or software packages:

	Data Processing was performe using code provided and described at:
	https://github.com/ArthurBlackburn/PtychoRunner

	This software in turn requires Matlab and has been tested to run successfully on Matlab versions 2019b and 2024b.
	
5. Was data derived from another source?
	A. no

6. Recommended citation for this dataset: 

	Arthur M. Blackburn, (2025), Electron Diffraction Data and Ptychographic Reconstructions Using a 20 keV Electron Beam on Supported Gold Nanoparticles, Federated Research Data Repository, doi: https://doi.org/10.20383/103.01357

---------------------
DATA & FILE OVERVIEW
---------------------

Root                                             
│  Readme.txt                                    Description of this dataset
│  DataAndRecipes.xls                            Spreadsheet file containing the metadata for the datasets and their reconstructions
│                                                
├─ Source                                        ** Source data for reconstructions **
│   │                                            
│   ├─ Experimental                              ** Directory containing source experimental data **
│   │   Au_on_aC_full_uncorr.mat                 Diffraction data from Au particles on amorphous carbon (aC), with no distortion correction, 90 by 90 set
│   │   Au_on_aC_subset_A_uncorr.mat             Diffraction data from Au particles on amorphous carbon (aC), with no distortion correction, sub-set A
│   │   Au_on_aC_subset_B_uncorr.mat             Diffraction data from Au particles on amorphous carbon (aC), with no distortion correction, sub-set B
│   │   Au_on_aC_full_cubic_corr1.mat            Diffraction data from Au particles on amorphous carbon (aC), after first round of distortion correction, 90 by 90 set
│   │   Au_on_aC_subset_A_cubic_corr1.mat        Diffraction data from Au particles on amorphous carbon (aC), after first round of distortion correction, sub-set A
│   │   Au_on_aC_subset_B_cubic_corr1.mat        Diffraction data from Au particles on amorphous carbon (aC), after first round of distortion correction, sub-set B
│   │   Au_on_aC_full_second_corr.mat            Diffraction data from Au particles on amorphous carbon (aC), after second (and final) of distortion correction, 90 by 90 set
│   │   Au_on_aC_subset_A_second_corr.mat        Diffraction data from Au particles on amorphous carbon (aC), after second (and final) of distortion correction, sub-set A
│   │   Au_on_aC_subset_B_second_corr.mat        Diffraction data from Au particles on amorphous carbon (aC), after second (and final) of distortion correction, sub-set B
│   │   Au_on_aC_full_posn.hdf5                  Positional data for diffraction data of Au particles on amorphous carbon (aC), 90 by 90 set
│   │   Au_on_aC_subset_A_posn.hdf5              Positional data for diffraction data of Au particles on amorphous carbon (aC), sub-set A
│   │   Au_on_aC_subset_B_posn.hdf5              Positional data for diffraction data of Au particles on amorphous carbon (aC), sub-set B
│   │   Au_on_aC_source.hp                       Source data for Au particles on amorphous carbon as collected from SU9000 with Azorus software, in hdf5 format
│   │   Au_on_MoS2_1024_corr.mat                 Diffraction data from Au islands on MoS2, with distortion correction
│   │   Au_on_MoS2_1024_corr_posn.hdf5           Positional data for diffraction data from Au particles on MoS2
│   │   MoS2_source.hp                           Source data for Au islands on MoS2 as collected from SU9000 with Azorus software, in hdf5 format
│   │   DistortionFit_cubic_corr1.mat            Fitting parameters extracted from fitting on first round of reconstructions (see Part_3_DistortionFit_and_Correct.m in Runner/TestScripts)
│   │   DistortionFit_cubic_corr2.mat            Fitting parameters extracted from fitting on second round of reconstructions (see Part_3_DistortionFit_and_Correct.m in Runner/TestScripts)
│   │         
│   └─ Synthetic                                                               ** Directory containing synthetic data **
│       │                                                                      
│       ├─ TestImage                                                           ** Directory containing data related to abstract test image **
│       │   test2048_complex_point_15105432_obj_padded_20mrad_diff.mat         Diffraction data for conditions similar to Au/MoS2 experiment, with 20 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_15105432_obj_padded_10mrad_diff.mat         Diffraction data for conditions similar to Au/MoS2 experiment, with 10 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_15105432_obj_padded_5mrad_diff.mat          Diffraction data for conditions similar to Au/MoS2 experiment, with 5 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_AuaC512n_obj_padded_15mrad_diff.mat         Diffraction data for conditions similar to Au - aC experiment, with 15 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_AuaC512n_obj_padded_10mrad_diff.mat         Diffraction data for conditions similar to Au - aC experiment, with 10 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_AuaC512n_obj_padded_8mrad_diff.mat          Diffraction data for conditions similar to Au - aC experiment, with 7.5 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_15105432_obj_padded_20mrad_positions.hdf5   Positional data for diffraction data for conditions similar to Au/MoS2 experiment, with 20 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_15105432_obj_padded_10mrad_positions.hdf5   Positional data for diffraction data for conditions similar to Au/MoS2 experiment, with 10 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_15105432_obj_padded_5mrad_positions.hdf5    Positional data for diffraction data for conditions similar to Au/MoS2 experiment, with 5 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_AuaC512n_obj_padded_15mrad_positions.hdf5   Positional data for diffraction data for conditions similar to Au - aC experiment, with 15 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_AuaC512n_obj_padded_10mrad_positions.hdf5   Positional data for diffraction data for conditions similar to Au - aC experiment, with 10 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_AuaC512n_obj_padded_8mrad_positions.hdf5    Positional data for diffraction data for conditions similar to Au - aC experiment, with 7.5 mrad semi-beam angle at 20 kV
│       │   Demo01_DiffractionData.mat                                         Diffraction data for a quick test to check on the reconstruction and software setup
│       │   Demo01_DiffractionData_Posns.hdf5                                  Poistional data for diffraction data for a quick test to check on the reconstruction and software setup
│       │   testobject_2048pix_padded.mat                                      The test object used to create the quick test and the other synthetic data above
│       │     
│       └─ AuMoS2                                       ** Directory containing data related to model of Au on 5 layer MoS2 **
│           README.md                                   Description of how to reproduce the data
│           generate.py                                 Script to produce the diffraction data using prismatique code from https://github.com/mrfitzpa/prismatique
│           probe-wavefunction.mat                      The wave function of the probe used in the model, at the model entry plane
│           hrtem-wavefunction-for-tilted-sample.mat    The exit wavefunction expected from a entry plane wave-function impinging on the tilted model
│           atomic-coords.xyz                           The x,y,z coordinates of the atoms in the model 
│           hrtem-wavefunction-for-no-tilt-sample.mat   The exit wavefunction expected from a entry plane wave-function impinging on the non-tilted model
│           probe-positions.hdf5                        Positional data for the modelled diffraction data
│           DPs-for-Cc-eq-0.mat                         The modelled diffraction data
│           DP1024-for-Cc-eq-0.mat                      The modelled diffraction data as above, but here zero-padded to be 1024-by-1024
│             
└─ Reconstructions                                      ** Directory containing ptychographic reconstructions produced from source data **
    AuaC_cubic_corr1_R_01_recipe002.mat                 Reconstruction from Au-aC data, full set, after first round of distortion correction, intermediate output (recipe002)  
    AuaC_cubic_corr1_R_01_recipe003.mat                 Reconstruction from Au-aC data, full set, after first round of distortion correction, intermediate output (recipe003)
    AuaC_cubic_corr1_R_01_recipe004.mat                 Reconstruction from Au-aC data, full set, after first round of distortion correction, intermediate output (recipe004)
    AuaC_cubic_corr1_R_01_recipe005.mat                 Reconstruction from Au-aC data, full set, after first round of distortion correction, intermediate output (recipe005)
    AuaC_cubic_corr1_R_01_recipe006.mat                 Reconstruction from Au-aC data, full set, after first round of distortion correction, final output (recipe006)
    AuaC_cubic_corr1_setA_R_01_recipe002.mat            Reconstruction from Au-aC data, set A, after first round of distortion correction, intermediate output (recipe002)
    AuaC_cubic_corr1_setA_R_01_recipe003.mat            Reconstruction from Au-aC data, set A, after first round of distortion correction, intermediate output (recipe003)
    AuaC_cubic_corr1_setA_R_01_recipe004.mat            Reconstruction from Au-aC data, set A, after first round of distortion correction, final output (recipe004)
    AuaC_cubic_corr1_setB_R_01_recipe002.mat            Reconstruction from Au-aC data, set B, after first round of distortion correction, intermediate output (recipe002)
    AuaC_cubic_corr1_setB_R_01_recipe003.mat            Reconstruction from Au-aC data, set B, after first round of distortion correction, intermediate output (recipe003)
    AuaC_cubic_corr1_setB_R_01_recipe004.mat            Reconstruction from Au-aC data, set B, after first round of distortion correction, final output (recipe004)
    AuaC_second_corr_R_01_recipe002.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, intermediate output (recipe002)
    AuaC_second_corr_R_01_recipe003.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, intermediate output (recipe003)
    AuaC_second_corr_R_01_recipe004.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, intermediate output (recipe004)
    AuaC_second_corr_R_01_recipe005.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, intermediate output (recipe005)
    AuaC_second_corr_R_01_recipe006.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, final output (recipe006)
    AuaC_second_corr_R_01_recipe017.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, intermediate output (recipe017)
    AuaC_second_corr_R_01_recipe018.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, intermediate output (recipe018)
    AuaC_second_corr_R_01_recipe019.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, intermediate output (recipe019)
    AuaC_second_corr_R_01_recipe020.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, intermediate output (recipe020)
    AuaC_second_corr_R_01_recipe021.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, final output (recipe021)
    AuaC_second_corr_R_01_recipe062.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe062)
    AuaC_second_corr_R_01_recipe063.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe063)
    AuaC_second_corr_R_01_recipe064.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe064)
    AuaC_second_corr_R_01_recipe065.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe065)
    AuaC_second_corr_R_01_recipe066.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Flat Init Phase, With FT Aperture, final output (recipe066)
    AuaC_second_corr_R_01_recipe068.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe068)
    AuaC_second_corr_R_01_recipe069.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe069)
    AuaC_second_corr_R_01_recipe070.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe070)
    AuaC_second_corr_R_01_recipe071.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe071)
    AuaC_second_corr_R_01_recipe072.mat                 Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, With FT Aperture, final output (recipe072)
    AuaC_second_corr_setA_R_01_recipe002.mat            Reconstruction from Au-aC data, set A, after second round of distortion correction, intermediate output (recipe002)
    AuaC_second_corr_setA_R_01_recipe003.mat            Reconstruction from Au-aC data, set A, after second round of distortion correction, intermediate output (recipe003)
    AuaC_second_corr_setA_R_01_recipe004.mat            Reconstruction from Au-aC data, set A, after second round of distortion correction, final output (recipe004)
    AuaC_second_corr_setA_R_01_recipe017.mat            Reconstruction from Au-aC data, set A, after second round of distortion correction, Random Init Phase, final output (recipe017)
    AuaC_second_corr_setA_R_01_recipe018.mat            Reconstruction from Au-aC data, set A, after second round of distortion correction, Random Init Phase, intermediate output (recipe018)
    AuaC_second_corr_setA_R_01_recipe019.mat            Reconstruction from Au-aC data, set A, after second round of distortion correction, Random Init Phase, intermediate output (recipe019)
    AuaC_second_corr_setA_R_01_recipe062.mat            Reconstruction from Au-aC data, set A, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe062)
    AuaC_second_corr_setA_R_01_recipe063.mat            Reconstruction from Au-aC data, set A, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe063)
    AuaC_second_corr_setA_R_01_recipe064.mat            Reconstruction from Au-aC data, set A, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe064)
    AuaC_second_corr_setA_R_01_recipe068.mat            Reconstruction from Au-aC data, set A, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe068)
    AuaC_second_corr_setA_R_01_recipe069.mat            Reconstruction from Au-aC data, set A, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe069)
    AuaC_second_corr_setA_R_01_recipe070.mat            Reconstruction from Au-aC data, set A, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe070)
    AuaC_second_corr_setB_R_01_recipe002.mat            Reconstruction from Au-aC data, set B, after second round of distortion correction, intermediate output (recipe002)
    AuaC_second_corr_setB_R_01_recipe003.mat            Reconstruction from Au-aC data, set B, after second round of distortion correction, intermediate output (recipe003)
    AuaC_second_corr_setB_R_01_recipe004.mat            Reconstruction from Au-aC data, set B, after second round of distortion correction, final output (recipe004)
    AuaC_second_corr_setB_R_01_recipe017.mat            Reconstruction from Au-aC data, set B, after second round of distortion correction, Random Init Phase, final output (recipe017)
    AuaC_second_corr_setB_R_01_recipe018.mat            Reconstruction from Au-aC data, set B, after second round of distortion correction, Random Init Phase, intermediate output (recipe018)
    AuaC_second_corr_setB_R_01_recipe019.mat            Reconstruction from Au-aC data, set B, after second round of distortion correction, Random Init Phase, intermediate output (recipe019)
    AuaC_second_corr_setB_R_01_recipe062.mat            Reconstruction from Au-aC data, set B, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe062)
    AuaC_second_corr_setB_R_01_recipe063.mat            Reconstruction from Au-aC data, set B, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe063)
    AuaC_second_corr_setB_R_01_recipe064.mat            Reconstruction from Au-aC data, set B, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe064)
    AuaC_second_corr_setB_R_01_recipe068.mat            Reconstruction from Au-aC data, set B, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe068)
    AuaC_second_corr_setB_R_01_recipe069.mat            Reconstruction from Au-aC data, set B, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe069)
    AuaC_second_corr_setB_R_01_recipe070.mat            Reconstruction from Au-aC data, set B, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe070)
    AuaC_uncorr_R_01_recipe002.mat                      Reconstruction from Au-aC data, full set, no distortion correction, intermediate output (recipe002)
    AuaC_uncorr_R_01_recipe003.mat                      Reconstruction from Au-aC data, full set, no distortion correction, intermediate output (recipe003)
    AuaC_uncorr_R_01_recipe004.mat                      Reconstruction from Au-aC data, full set, no distortion correction, intermediate output (recipe004)
    AuaC_uncorr_R_01_recipe005.mat                      Reconstruction from Au-aC data, full set, no distortion correction, intermediate output (recipe005)
    AuaC_uncorr_R_01_recipe006.mat                      Reconstruction from Au-aC data, full set, no distortion correction, final output (recipe006)
    AuaC_uncorr_setA_R_01_recipe002.mat                 Reconstruction from Au-aC data, set A, no distortion correction, intermediate output (recipe002)
    AuaC_uncorr_setA_R_01_recipe003.mat                 Reconstruction from Au-aC data, set A, no distortion correction, intermediate output (recipe003)
    AuaC_uncorr_setA_R_01_recipe004.mat                 Reconstruction from Au-aC data, set A, no distortion correction, final output (recipe004)
    AuaC_uncorr_setB_R_01_recipe002.mat                 Reconstruction from Au-aC data, set B, no distortion correction, intermediate output (recipe002)
    AuaC_uncorr_setB_R_01_recipe003.mat                 Reconstruction from Au-aC data, set B, no distortion correction, intermediate output (recipe003)
    AuaC_uncorr_setB_R_01_recipe004.mat                 Reconstruction from Au-aC data, set B, no distortion correction, final output (recipe004)
    MoS2_Au_model_R_01_recipe074.mat                    Reconstruction from Au - MoS2 model synthetic data, intermediate output (recipe074)
    MoS2_Au_model_R_01_recipe075.mat                    Reconstruction from Au - MoS2 model synthetic data, intermediate output (recipe075)
    MoS2_Au_model_R_01_recipe076.mat                    Reconstruction from Au - MoS2 model synthetic data, intermediate output (recipe076)
    MoS2_Au_model_R_01_recipe077.mat                    Reconstruction from Au - MoS2 model synthetic data, intermediate output (recipe077)
    MoS2_Au_model_R_01_recipe078.mat                    Reconstruction from Au - MoS2 model synthetic data, intermediate output (recipe078)
    MoS2_Au_model_R_01_recipe079.mat                    Reconstruction from Au - MoS2 model synthetic data, intermediate output (recipe079)
    MoS2_corrected_R_01_recipe008.mat                   Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe008)
    MoS2_corrected_R_01_recipe009.mat                   Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe009)
    MoS2_corrected_R_01_recipe010.mat                   Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe010)
    MoS2_corrected_R_01_recipe011.mat                   Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe011)
    MoS2_corrected_R_01_recipe012.mat                   Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe012)
    MoS2_corrected_R_01_recipe013.mat                   Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe013)
    MoS2_corrected_R_01_recipe014.mat                   Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe014)
    MoS2_corrected_R_01_recipe015.mat                   Reconstruction from Au - MoS2 experimental data, distortion correction applied, final output (recipe015)
    synth_aCconds_5m_R_01_recipe056.mat                 Reconstruction from synethetic single slice, conditions similar to Au-aC though with 5 mrad beam, intermediate output (recipe056)
    synth_aCconds_5m_R_01_recipe057.mat                 Reconstruction from synethetic single slice, conditions similar to Au-aC though with 5 mrad beam, intermediate output (recipe057)
    synth_aCconds_5m_R_01_recipe058.mat                 Reconstruction from synethetic single slice, conditions similar to Au-aC though with 5 mrad beam, intermediate output (recipe058)
    synth_aCconds_5m_R_01_recipe059.mat                 Reconstruction from synethetic single slice, conditions similar to Au-aC though with 5 mrad beam, intermediate output (recipe059)
    synth_aCconds_5m_R_01_recipe060.mat                 Reconstruction from synethetic single slice, conditions similar to Au-aC though with 5 mrad beam, final output (recipe060)
    synth_aCconds_10m_R_01_recipe050.mat                Reconstruction from synethetic single slice, conditions similar to Au-aC though with 10 mrad beam, intermediate output (recipe050)
    synth_aCconds_10m_R_01_recipe051.mat                Reconstruction from synethetic single slice, conditions similar to Au-aC though with 10 mrad beam, intermediate output (recipe051)
    synth_aCconds_10m_R_01_recipe052.mat                Reconstruction from synethetic single slice, conditions similar to Au-aC though with 10 mrad beam, intermediate output (recipe052)
    synth_aCconds_10m_R_01_recipe053.mat                Reconstruction from synethetic single slice, conditions similar to Au-aC though with 10 mrad beam, intermediate output (recipe053)
    synth_aCconds_10m_R_01_recipe054.mat                Reconstruction from synethetic single slice, conditions similar to Au-aC though with 10 mrad beam, final output (recipe054)
    synth_aCconds_15m_R_01_recipe044.mat                Reconstruction from synethetic single slice, conditions similar to Au-aC though with 15 mrad beam, intermediate output (recipe044)
    synth_aCconds_15m_R_01_recipe045.mat                Reconstruction from synethetic single slice, conditions similar to Au-aC though with 15 mrad beam, intermediate output (recipe045)
    synth_aCconds_15m_R_01_recipe046.mat                Reconstruction from synethetic single slice, conditions similar to Au-aC though with 15 mrad beam, intermediate output (recipe046)
    synth_aCconds_15m_R_01_recipe047.mat                Reconstruction from synethetic single slice, conditions similar to Au-aC though with 15 mrad beam, intermediate output (recipe047)
    synth_aCconds_15m_R_01_recipe048.mat                Reconstruction from synethetic single slice, conditions similar to Au-aC though with 15 mrad beam, final output (recipe048)
    synth_demo01_R_01_recipe023.mat                     Reconstruction from simple test data
    synth_MoS2conds_5m_R_01_recipe038.mat               Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 7.5 mrad beam, intermediate output (recipe038)
    synth_MoS2conds_5m_R_01_recipe039.mat               Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 7.5 mrad beam, intermediate output (recipe039)
    synth_MoS2conds_5m_R_01_recipe040.mat               Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 7.5 mrad beam, intermediate output (recipe040)
    synth_MoS2conds_5m_R_01_recipe041.mat               Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 7.5 mrad beam, intermediate output (recipe041)
    synth_MoS2conds_5m_R_01_recipe042.mat               Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 7.5 mrad beam, final output (recipe042)
    synth_MoS2conds_10m_R_01_recipe032.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 10 mrad beam, intermediate output (recipe032)
    synth_MoS2conds_10m_R_01_recipe033.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 10 mrad beam, intermediate output (recipe033)
    synth_MoS2conds_10m_R_01_recipe034.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 10 mrad beam, intermediate output (recipe034)
    synth_MoS2conds_10m_R_01_recipe035.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 10 mrad beam, intermediate output (recipe035)
    synth_MoS2conds_10m_R_01_recipe036.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 10 mrad beam, final output (recipe036)
    synth_MoS2conds_20m_R_01_recipe026.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 15 mrad beam, intermediate output (recipe026)
    synth_MoS2conds_20m_R_01_recipe027.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 15 mrad beam, intermediate output (recipe027)
    synth_MoS2conds_20m_R_01_recipe028.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 15 mrad beam, intermediate output (recipe028)
    synth_MoS2conds_20m_R_01_recipe029.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 15 mrad beam, intermediate output (recipe029)
    synth_MoS2conds_20m_R_01_recipe030.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 15 mrad beam, final output (recipe030)
    AuaC_cubic_corr1_R_01_recipe002_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after first round of distortion correction, intermediate output (recipe002)
    AuaC_cubic_corr1_R_01_recipe003_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after first round of distortion correction, intermediate output (recipe003)
    AuaC_cubic_corr1_R_01_recipe004_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after first round of distortion correction, intermediate output (recipe004)
    AuaC_cubic_corr1_R_01_recipe005_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after first round of distortion correction, intermediate output (recipe005)
    AuaC_cubic_corr1_R_01_recipe006_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after first round of distortion correction, final output (recipe006)
    AuaC_cubic_corr1_setA_R_01_recipe002_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after first round of distortion correction, intermediate output (recipe002)  
    AuaC_cubic_corr1_setA_R_01_recipe003_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after first round of distortion correction, intermediate output (recipe003)  
    AuaC_cubic_corr1_setA_R_01_recipe004_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after first round of distortion correction, final output (recipe004)  
    AuaC_cubic_corr1_setB_R_01_recipe002_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after first round of distortion correction, intermediate output (recipe002)  
    AuaC_cubic_corr1_setB_R_01_recipe003_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after first round of distortion correction, intermediate output (recipe003)  
    AuaC_cubic_corr1_setB_R_01_recipe004_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after first round of distortion correction, final output (recipe004)  
    AuaC_second_corr_R_01_recipe002_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, intermediate output (recipe002)
    AuaC_second_corr_R_01_recipe003_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, intermediate output (recipe003)
    AuaC_second_corr_R_01_recipe004_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, intermediate output (recipe004)
    AuaC_second_corr_R_01_recipe005_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, intermediate output (recipe005)
    AuaC_second_corr_R_01_recipe006_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, final output (recipe006)
    AuaC_second_corr_R_01_recipe017_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, intermediate output (recipe017)
    AuaC_second_corr_R_01_recipe018_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, intermediate output (recipe018)
    AuaC_second_corr_R_01_recipe019_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, intermediate output (recipe019)
    AuaC_second_corr_R_01_recipe020_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, intermediate output (recipe020)
    AuaC_second_corr_R_01_recipe021_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, final output (recipe021)
    AuaC_second_corr_R_01_recipe062_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe062)
    AuaC_second_corr_R_01_recipe063_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe063)
    AuaC_second_corr_R_01_recipe064_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe064)
    AuaC_second_corr_R_01_recipe065_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe065)
    AuaC_second_corr_R_01_recipe066_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Flat Init Phase, With FT Aperture, final output (recipe066)
    AuaC_second_corr_R_01_recipe068_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe068)
    AuaC_second_corr_R_01_recipe069_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe069)
    AuaC_second_corr_R_01_recipe070_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe070)
    AuaC_second_corr_R_01_recipe071_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe071)
    AuaC_second_corr_R_01_recipe072_summary.tif         Summary image of Reconstruction from Au-aC data, full set, after second round of distortion correction, Random Init Phase, With FT Aperture, final output (recipe072)
    AuaC_second_corr_setA_R_01_recipe002_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after second round of distortion correction, intermediate output (recipe002)  
    AuaC_second_corr_setA_R_01_recipe003_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after second round of distortion correction, intermediate output (recipe003)  
    AuaC_second_corr_setA_R_01_recipe004_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after second round of distortion correction, final output (recipe004)  
    AuaC_second_corr_setA_R_01_recipe017_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after second round of distortion correction, Random Init Phase, final output (recipe017)  
    AuaC_second_corr_setA_R_01_recipe018_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after second round of distortion correction, Random Init Phase, intermediate output (recipe018)  
    AuaC_second_corr_setA_R_01_recipe019_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after second round of distortion correction, Random Init Phase, intermediate output (recipe019)  
    AuaC_second_corr_setA_R_01_recipe062_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe062)  
    AuaC_second_corr_setA_R_01_recipe063_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe063)  
    AuaC_second_corr_setA_R_01_recipe064_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe064)  
    AuaC_second_corr_setA_R_01_recipe068_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe068)  
    AuaC_second_corr_setA_R_01_recipe069_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe069)  
    AuaC_second_corr_setA_R_01_recipe070_summary.tif    Summary image of Reconstruction from Au-aC data, set A, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe070)  
    AuaC_second_corr_setB_R_01_recipe002_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after second round of distortion correction, intermediate output (recipe002)  
    AuaC_second_corr_setB_R_01_recipe003_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after second round of distortion correction, intermediate output (recipe003)  
    AuaC_second_corr_setB_R_01_recipe004_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after second round of distortion correction, final output (recipe004)  
    AuaC_second_corr_setB_R_01_recipe017_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after second round of distortion correction, Random Init Phase, final output (recipe017)  
    AuaC_second_corr_setB_R_01_recipe018_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after second round of distortion correction, Random Init Phase, intermediate output (recipe018)  
    AuaC_second_corr_setB_R_01_recipe019_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after second round of distortion correction, Random Init Phase, intermediate output (recipe019)  
    AuaC_second_corr_setB_R_01_recipe062_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe062)  
    AuaC_second_corr_setB_R_01_recipe063_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe063)  
    AuaC_second_corr_setB_R_01_recipe064_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after second round of distortion correction, Flat Init Phase, With FT Aperture, intermediate output (recipe064)  
    AuaC_second_corr_setB_R_01_recipe068_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe068)  
    AuaC_second_corr_setB_R_01_recipe069_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe069)  
    AuaC_second_corr_setB_R_01_recipe070_summary.tif    Summary image of Reconstruction from Au-aC data, set B, after second round of distortion correction, Random Init Phase, With FT Aperture, intermediate output (recipe070)  
    AuaC_uncorr_R_01_recipe002_summary.tif              Summary image of Reconstruction from Au-aC data, full set, no distortion correction, intermediate output (recipe002)
    AuaC_uncorr_R_01_recipe003_summary.tif              Summary image of Reconstruction from Au-aC data, full set, no distortion correction, intermediate output (recipe003)
    AuaC_uncorr_R_01_recipe004_summary.tif              Summary image of Reconstruction from Au-aC data, full set, no distortion correction, intermediate output (recipe004)
    AuaC_uncorr_R_01_recipe005_summary.tif              Summary image of Reconstruction from Au-aC data, full set, no distortion correction, intermediate output (recipe005)
    AuaC_uncorr_R_01_recipe006_summary.tif              Summary image of Reconstruction from Au-aC data, full set, no distortion correction, final output (recipe006)
    AuaC_uncorr_setA_R_01_recipe002_summary.tif         Summary image of Reconstruction from Au-aC data, set A, no distortion correction, intermediate output (recipe002)
    AuaC_uncorr_setA_R_01_recipe003_summary.tif         Summary image of Reconstruction from Au-aC data, set A, no distortion correction, intermediate output (recipe003)
    AuaC_uncorr_setA_R_01_recipe004_summary.tif         Summary image of Reconstruction from Au-aC data, set A, no distortion correction, final output (recipe004)
    AuaC_uncorr_setB_R_01_recipe002_summary.tif         Summary image of Reconstruction from Au-aC data, set B, no distortion correction, intermediate output (recipe002)
    AuaC_uncorr_setB_R_01_recipe003_summary.tif         Summary image of Reconstruction from Au-aC data, set B, no distortion correction, intermediate output (recipe003)
    AuaC_uncorr_setB_R_01_recipe004_summary.tif         Summary image of Reconstruction from Au-aC data, set B, no distortion correction, final output (recipe004)
    MoS2_Au_model_R_01_recipe074_summary.tif            Summary image of Reconstruction from Au - MoS2 model synthetic data, intermediate output (recipe074)
    MoS2_Au_model_R_01_recipe075_summary.tif            Summary image of Reconstruction from Au - MoS2 model synthetic data, intermediate output (recipe075)
    MoS2_Au_model_R_01_recipe076_summary.tif            Summary image of Reconstruction from Au - MoS2 model synthetic data, intermediate output (recipe076)
    MoS2_Au_model_R_01_recipe077_summary.tif            Summary image of Reconstruction from Au - MoS2 model synthetic data, intermediate output (recipe077)
    MoS2_Au_model_R_01_recipe078_summary.tif            Summary image of Reconstruction from Au - MoS2 model synthetic data, intermediate output (recipe078)
    MoS2_Au_model_R_01_recipe079_summary.tif            Summary image of Reconstruction from Au - MoS2 model synthetic data, intermediate output (recipe079)
    MoS2_corrected_R_01_recipe008_summary.tif           Summary image of Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe008)
    MoS2_corrected_R_01_recipe009_summary.tif           Summary image of Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe009)
    MoS2_corrected_R_01_recipe010_summary.tif           Summary image of Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe010)
    MoS2_corrected_R_01_recipe011_summary.tif           Summary image of Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe011)
    MoS2_corrected_R_01_recipe012_summary.tif           Summary image of Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe012)
    MoS2_corrected_R_01_recipe013_summary.tif           Summary image of Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe013)
    MoS2_corrected_R_01_recipe014_summary.tif           Summary image of Reconstruction from Au - MoS2 experimental data, distortion correction applied, intermediate output (recipe014)
    MoS2_corrected_R_01_recipe015_summary.tif           Summary image of Reconstruction from Au - MoS2 experimental data, distortion correction applied, final output (recipe015)
    synth_aCconds_5m_R_01_recipe056_summary.tif         Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 5 mrad beam, intermediate output (recipe056)
    synth_aCconds_5m_R_01_recipe057_summary.tif         Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 5 mrad beam, intermediate output (recipe057)
    synth_aCconds_5m_R_01_recipe058_summary.tif         Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 5 mrad beam, intermediate output (recipe058)
    synth_aCconds_5m_R_01_recipe059_summary.tif         Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 5 mrad beam, intermediate output (recipe059)
    synth_aCconds_5m_R_01_recipe060_summary.tif         Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 5 mrad beam, final output (recipe060)
    synth_aCconds_10m_R_01_recipe050_summary.tif        Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 10 mrad beam, intermediate output (recipe050)
    synth_aCconds_10m_R_01_recipe051_summary.tif        Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 10 mrad beam, intermediate output (recipe051)
    synth_aCconds_10m_R_01_recipe052_summary.tif        Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 10 mrad beam, intermediate output (recipe052)
    synth_aCconds_10m_R_01_recipe053_summary.tif        Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 10 mrad beam, intermediate output (recipe053)
    synth_aCconds_10m_R_01_recipe054_summary.tif        Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 10 mrad beam, final output (recipe054)
    synth_aCconds_15m_R_01_recipe044_summary.tif        Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 15 mrad beam, intermediate output (recipe044)
    synth_aCconds_15m_R_01_recipe045_summary.tif        Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 15 mrad beam, intermediate output (recipe045)
    synth_aCconds_15m_R_01_recipe046_summary.tif        Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 15 mrad beam, intermediate output (recipe046)
    synth_aCconds_15m_R_01_recipe047_summary.tif        Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 15 mrad beam, intermediate output (recipe047)
    synth_aCconds_15m_R_01_recipe048_summary.tif        Summary image of Reconstruction from synethetic single slice, conditions similar to Au-aC though with 15 mrad beam, final output (recipe048)
    synth_demo01_R_01_recipe023_summary.tif             Summary image of Reconstruction from simple test data
    synth_MoS2conds_5m_R_01_recipe038_summary.tif       Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 7.5 mrad beam, intermediate output (recipe038)
    synth_MoS2conds_5m_R_01_recipe039_summary.tif       Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 7.5 mrad beam, intermediate output (recipe039)
    synth_MoS2conds_5m_R_01_recipe040_summary.tif       Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 7.5 mrad beam, intermediate output (recipe040)
    synth_MoS2conds_5m_R_01_recipe041_summary.tif       Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 7.5 mrad beam, intermediate output (recipe041)
    synth_MoS2conds_5m_R_01_recipe042_summary.tif       Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 7.5 mrad beam, final output (recipe042)
    synth_MoS2conds_10m_R_01_recipe032_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 10 mrad beam, intermediate output (recipe032)
    synth_MoS2conds_10m_R_01_recipe033_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 10 mrad beam, intermediate output (recipe033)
    synth_MoS2conds_10m_R_01_recipe034_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 10 mrad beam, intermediate output (recipe034)
    synth_MoS2conds_10m_R_01_recipe035_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 10 mrad beam, intermediate output (recipe035)
    synth_MoS2conds_10m_R_01_recipe036_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 10 mrad beam, final output (recipe036)
    synth_MoS2conds_20m_R_01_recipe026_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 15 mrad beam, intermediate output (recipe026)
    synth_MoS2conds_20m_R_01_recipe027_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 15 mrad beam, intermediate output (recipe027)
    synth_MoS2conds_20m_R_01_recipe028_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 15 mrad beam, intermediate output (recipe028)
    synth_MoS2conds_20m_R_01_recipe029_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 15 mrad beam, intermediate output (recipe029)
    synth_MoS2conds_20m_R_01_recipe030_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS2 though with 15 mrad beam, final output (recipe030)                              


---------------------------
METHODOLOGICAL INFORMATION
---------------------------

1. Description of methods used for collection/generation of data:

The methods of data collection are described in the associated publication, linked to at https://github.com/ArthurBlackburn/PtychoRunner. 

2. Methods for processing the data: 

The processing of the data is also described in the associated publication, and within the associate code provided at https://github.com/ArthurBlackburn/PtychoRunner. 

3. Instrument- or software-specific information needed to interpret the data:

Experimental data was collected on a modified Hitachi SU9000 electron microscope, fitted with Dectris Quadro detector optimized for low energy electron detection (see Tinti, G. et al. The EIGER detector for low-energy electron microscopy and photoemission electron microscopy. Journal of Synchrotron Radiation 24, 963-974 (2017). https://doi.org/10.1107/S1600577517009109). The *.mat files above store data in Matlab v7.3 format. Version 7.3 MAT-files use an HDF5-based format that stores data in compressed chunks. The *.hp and %.hdf5 files above are HDF5 format, and can be opened using standard HDF5 libraries, such as described at https://github.com/HDFGroup/.

4. Standards and calibration information, if appropriate: 

The distortion corrected gold on amorhphous diffraction datasets have been calirbrated and corrected against the known lattice information for gold.


