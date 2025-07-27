This README.txt file was generated on 2025-07-18 by Arthur M. Blackburn

--------------------
GENERAL INFORMATION
--------------------

1. Title of Dataset: Electron Diffraction Data and Ptychographic Reconstructions Using a 20 keV Electron Beam on Supported Gold Nanoparticles

2. Author Information
	A. Principal Investigator Contact Information
		Name: Arthur M. BLACKBURN
		Institution: University of Victoria
		Email: ablackbu@uvic.ca
		Alternative Email: arthur.blackburn@gmail.com
	
	B. Co-Authors
		Name: Cristina Cordoba
		Institution: University of Victoria  
		
		Name: Matthew R. Fitzpatrick
		Institution: University of Victoria
		
		Name: Robert A. McLeod
		Institution: Hitachi High-Tech Canada Inc.

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
	A. No

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
│   │   Au_on_MoS₂_1024_corr.mat                 Diffraction data from Au islands on MoS₂, with distortion correction
│   │   Au_on_MoS₂_1024_corr_posn.hdf5           Positional data for diffraction data from Au particles on MoS₂
│   │   MoS₂_source.hp                           Source data for Au islands on MoS₂ as collected from SU9000 with Azorus software, in hdf5 format
│   │   DistortionFit_cubic_corr1.mat            Fitting parameters extracted from fitting on first round of reconstructions (see Part_3_DistortionFit_and_Correct.m in Runner/TestScripts)
│   │   DistortionFit_cubic_corr2.mat            Fitting parameters extracted from fitting on second round of reconstructions (see Part_3_DistortionFit_and_Correct.m in Runner/TestScripts)
│   │         
│   └─ Synthetic                                                               ** Directory containing synthetic data **
│       │                                                                      
│       ├─ TestImage                                                           ** Directory containing data related to abstract test image **
│       │   test2048_complex_point_15105432_obj_padded_20mrad_diff.mat         Diffraction data for conditions similar to Au/MoS₂ experiment, with 20 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_15105432_obj_padded_10mrad_diff.mat         Diffraction data for conditions similar to Au/MoS₂ experiment, with 10 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_15105432_obj_padded_5mrad_diff.mat          Diffraction data for conditions similar to Au/MoS₂ experiment, with 5 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_AuaC512n_obj_padded_15mrad_diff.mat         Diffraction data for conditions similar to Au - aC experiment, with 15 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_AuaC512n_obj_padded_10mrad_diff.mat         Diffraction data for conditions similar to Au - aC experiment, with 10 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_AuaC512n_obj_padded_8mrad_diff.mat          Diffraction data for conditions similar to Au - aC experiment, with 7.5 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_15105432_obj_padded_20mrad_positions.hdf5   Positional data for diffraction data for conditions similar to Au/MoS₂ experiment, with 20 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_15105432_obj_padded_10mrad_positions.hdf5   Positional data for diffraction data for conditions similar to Au/MoS₂ experiment, with 10 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_15105432_obj_padded_5mrad_positions.hdf5    Positional data for diffraction data for conditions similar to Au/MoS₂ experiment, with 5 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_AuaC512n_obj_padded_15mrad_positions.hdf5   Positional data for diffraction data for conditions similar to Au - aC experiment, with 15 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_AuaC512n_obj_padded_10mrad_positions.hdf5   Positional data for diffraction data for conditions similar to Au - aC experiment, with 10 mrad semi-beam angle at 20 kV
│       │   test2048_complex_point_AuaC512n_obj_padded_8mrad_positions.hdf5    Positional data for diffraction data for conditions similar to Au - aC experiment, with 7.5 mrad semi-beam angle at 20 kV
│       │   Demo01_DiffractionData.mat                                         Diffraction data for a quick test to check on the reconstruction and software setup
│       │   Demo01_DiffractionData_Posns.hdf5                                  Poistional data for diffraction data for a quick test to check on the reconstruction and software setup
│       │   testobject_2048pix_padded.mat                                      The test object used to create the quick test and the other synthetic data above
│       │     
│       └─ AuMoS₂                                       ** Directory containing data related to model of Au on 5 layer MoS₂ **
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
    MoS₂_Au_model_R_01_recipe074.mat                    Reconstruction from Au - MoS₂ model synthetic data, intermediate output (recipe074)
    MoS₂_Au_model_R_01_recipe075.mat                    Reconstruction from Au - MoS₂ model synthetic data, intermediate output (recipe075)
    MoS₂_Au_model_R_01_recipe076.mat                    Reconstruction from Au - MoS₂ model synthetic data, intermediate output (recipe076)
    MoS₂_Au_model_R_01_recipe077.mat                    Reconstruction from Au - MoS₂ model synthetic data, intermediate output (recipe077)
    MoS₂_Au_model_R_01_recipe078.mat                    Reconstruction from Au - MoS₂ model synthetic data, intermediate output (recipe078)
    MoS₂_Au_model_R_01_recipe079.mat                    Reconstruction from Au - MoS₂ model synthetic data, intermediate output (recipe079)
    MoS₂_corrected_R_01_recipe008.mat                   Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe008)
    MoS₂_corrected_R_01_recipe009.mat                   Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe009)
    MoS₂_corrected_R_01_recipe010.mat                   Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe010)
    MoS₂_corrected_R_01_recipe011.mat                   Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe011)
    MoS₂_corrected_R_01_recipe012.mat                   Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe012)
    MoS₂_corrected_R_01_recipe013.mat                   Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe013)
    MoS₂_corrected_R_01_recipe014.mat                   Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe014)
    MoS₂_corrected_R_01_recipe015.mat                   Reconstruction from Au - MoS₂ experimental data, distortion correction applied, final output (recipe015)
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
    synth_MoS₂conds_5m_R_01_recipe038.mat               Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 7.5 mrad beam, intermediate output (recipe038)
    synth_MoS₂conds_5m_R_01_recipe039.mat               Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 7.5 mrad beam, intermediate output (recipe039)
    synth_MoS₂conds_5m_R_01_recipe040.mat               Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 7.5 mrad beam, intermediate output (recipe040)
    synth_MoS₂conds_5m_R_01_recipe041.mat               Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 7.5 mrad beam, intermediate output (recipe041)
    synth_MoS₂conds_5m_R_01_recipe042.mat               Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 7.5 mrad beam, final output (recipe042)
    synth_MoS₂conds_10m_R_01_recipe032.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 10 mrad beam, intermediate output (recipe032)
    synth_MoS₂conds_10m_R_01_recipe033.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 10 mrad beam, intermediate output (recipe033)
    synth_MoS₂conds_10m_R_01_recipe034.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 10 mrad beam, intermediate output (recipe034)
    synth_MoS₂conds_10m_R_01_recipe035.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 10 mrad beam, intermediate output (recipe035)
    synth_MoS₂conds_10m_R_01_recipe036.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 10 mrad beam, final output (recipe036)
    synth_MoS₂conds_20m_R_01_recipe026.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 15 mrad beam, intermediate output (recipe026)
    synth_MoS₂conds_20m_R_01_recipe027.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 15 mrad beam, intermediate output (recipe027)
    synth_MoS₂conds_20m_R_01_recipe028.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 15 mrad beam, intermediate output (recipe028)
    synth_MoS₂conds_20m_R_01_recipe029.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 15 mrad beam, intermediate output (recipe029)
    synth_MoS₂conds_20m_R_01_recipe030.mat              Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 15 mrad beam, final output (recipe030)
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
    MoS₂_Au_model_R_01_recipe074_summary.tif            Summary image of Reconstruction from Au - MoS₂ model synthetic data, intermediate output (recipe074)
    MoS₂_Au_model_R_01_recipe075_summary.tif            Summary image of Reconstruction from Au - MoS₂ model synthetic data, intermediate output (recipe075)
    MoS₂_Au_model_R_01_recipe076_summary.tif            Summary image of Reconstruction from Au - MoS₂ model synthetic data, intermediate output (recipe076)
    MoS₂_Au_model_R_01_recipe077_summary.tif            Summary image of Reconstruction from Au - MoS₂ model synthetic data, intermediate output (recipe077)
    MoS₂_Au_model_R_01_recipe078_summary.tif            Summary image of Reconstruction from Au - MoS₂ model synthetic data, intermediate output (recipe078)
    MoS₂_Au_model_R_01_recipe079_summary.tif            Summary image of Reconstruction from Au - MoS₂ model synthetic data, intermediate output (recipe079)
    MoS₂_corrected_R_01_recipe008_summary.tif           Summary image of Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe008)
    MoS₂_corrected_R_01_recipe009_summary.tif           Summary image of Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe009)
    MoS₂_corrected_R_01_recipe010_summary.tif           Summary image of Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe010)
    MoS₂_corrected_R_01_recipe011_summary.tif           Summary image of Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe011)
    MoS₂_corrected_R_01_recipe012_summary.tif           Summary image of Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe012)
    MoS₂_corrected_R_01_recipe013_summary.tif           Summary image of Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe013)
    MoS₂_corrected_R_01_recipe014_summary.tif           Summary image of Reconstruction from Au - MoS₂ experimental data, distortion correction applied, intermediate output (recipe014)
    MoS₂_corrected_R_01_recipe015_summary.tif           Summary image of Reconstruction from Au - MoS₂ experimental data, distortion correction applied, final output (recipe015)
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
    synth_MoS₂conds_5m_R_01_recipe038_summary.tif       Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 7.5 mrad beam, intermediate output (recipe038)
    synth_MoS₂conds_5m_R_01_recipe039_summary.tif       Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 7.5 mrad beam, intermediate output (recipe039)
    synth_MoS₂conds_5m_R_01_recipe040_summary.tif       Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 7.5 mrad beam, intermediate output (recipe040)
    synth_MoS₂conds_5m_R_01_recipe041_summary.tif       Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 7.5 mrad beam, intermediate output (recipe041)
    synth_MoS₂conds_5m_R_01_recipe042_summary.tif       Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 7.5 mrad beam, final output (recipe042)
    synth_MoS₂conds_10m_R_01_recipe032_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 10 mrad beam, intermediate output (recipe032)
    synth_MoS₂conds_10m_R_01_recipe033_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 10 mrad beam, intermediate output (recipe033)
    synth_MoS₂conds_10m_R_01_recipe034_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 10 mrad beam, intermediate output (recipe034)
    synth_MoS₂conds_10m_R_01_recipe035_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 10 mrad beam, intermediate output (recipe035)
    synth_MoS₂conds_10m_R_01_recipe036_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 10 mrad beam, final output (recipe036)
    synth_MoS₂conds_20m_R_01_recipe026_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 15 mrad beam, intermediate output (recipe026)
    synth_MoS₂conds_20m_R_01_recipe027_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 15 mrad beam, intermediate output (recipe027)
    synth_MoS₂conds_20m_R_01_recipe028_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 15 mrad beam, intermediate output (recipe028)
    synth_MoS₂conds_20m_R_01_recipe029_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 15 mrad beam, intermediate output (recipe029)
    synth_MoS₂conds_20m_R_01_recipe030_summary.tif      Summary image of Reconstruction from synethetic single slice, conditions similar to Au-MoS₂ though with 15 mrad beam, final output (recipe030)                              


---------------------------
METHODOLOGICAL INFORMATION
---------------------------

1. Description of methods used for collection/generation of data:

For the Au/aC experiments we used a combined gold and graphitized carbon on amorphous carbon sample from Ted Pella Inc. (CA, USA), part number #638. This sample underwent approximately 2 minutes of cleaning in a Hitachi High-Tech ZoneTEM ozone – deep UV cleaning system at its lowest pressure, to remove residual organic contamination and to slightly thin the aC support film by ~ 2 nm. This gave an underlying aC thickness of ~15 nm, as measured by electron energy loss spectroscopy in a separate electron microscope at 200 keV, by analyzing the ratio of zero-loss peak intensity to the total intensity. The preparation process for the Au/MoS₂ sample is described elsewhere[1]. However, here we determined that the MoS₂ was composed of 5 mono-layers which implies a thickness of 3.1 nm, using model-fitting between simulated and experimental convergent beam electron diffraction patterns, obtained by 4D-STEM [2]. The gold islands on the MoS₂ have diameters ranging from about 5 – 20 nm, with thickness of a smaller island estimated to be approximately 4 nm through analysis of the phase of our ptychographic reconstruction and matching to simulation.

Experiments were performed on a modified Hitachi SU9000 scanning electron microscope (SEM). This instrument can produce a maximum beam energy of 30 keV, though here it was used at 20 keV. This SEM is also equipped with a cold field-emission electron gun (CFEG) stabilized with a non-evaporative getter pump. Samples are placed in the immersion-objective magnetic-lens gap using a conventional side-entry style transmission electron microscope (TEM) holder. For our work we used a single tilt TEM sample holder. Compared to a double-tilt holder, this made achieving an on-axis condition impractical with the MoS₂ sample, particularly given that we aimed to minimize its time under electron beam irradiation to reduce possible beam damage and minimize beam-induced build up of surface contamination. The instrument was modified from a standard design to include a small projector lens beneath the objective lens. This additional lens improves control of the camera length (magnification of the diffraction plane) which is realized in combination with varying the objective lens current and the sample z-position (height) within the objective lens gap.
 
Beneath the retractable bright field detector of the SU9000, we added a hybrid type pixel array direct electron detector. This detector is based on the EIGER detector design[3], as provided in a Quadro family camera from Dectris AG, Switzerland. For this work the detector was customized for low energy electron (6 - 30keV) performance by removing the uppermost aluminium layer, thus eliminating a dead layer in the detector which becomes increasingly significant at low energies. The detector collects 512-by-512 pixel-count images (where each pixel is 75 μm square) and can operate at a full-frame rate of 4500 Hz (2250 Hz) when operating at 8-bit (or 16-bit) pixel depth. However, for our experiments we used a frame rate of 200 Hz, with no deadtime between frames, giving a 5 ms exposure for each collected diffraction pattern. At our chosen counting threshold, we determined a mean gain of 0.85 counts/e-.
 
For the gold on amorphous carbon (Au/aC) data acquisition our mean beam step size was 0.5 nm based on a square grid of 90 by 90 points, whereas for the Au/MoS₂ sample a step size of 0.9 nm was used on a 45 by 45 grid. This relates to 45-by-45 nm and 40.5-by-40.5 nm fields of view, acquired in 40.5 and 10 seconds respectively for the Au/aC and Au/MoS₂ datasets. Pseudorandom, but known, offsets from this grid no greater than ± 50% of the step size were added to the otherwise regular grid. The mean effective overlap of the illumination regions determined from the probe reconstructions were thus determined to be 0.875 and 0.74, given the effective illumination diameters from the probe reconstructions of approximately 4 nm and 3.5 nm for Au/aC and Au/MoS₂, respectively.
 
The electron beam current, measured using a pre-specimen insertable Faraday cup, was 85 pA for the Au/aC and 40 pA Au/MoS₂ experiments. Thus, given the fields of view, the mean electron dose received by these samples were thus 1.06 × 10^5 e/Å^2 and 1.56 × 10^4 e/Å^2, respectively. For the Au/aC data, the mean pixel count received on the Quadro detector was 5.035 counts/pix. For our 5 ms frame time, this gives an average rate of 2.64 × 10^8 counts/sec in total across the detector. With our detector giving 0.85 counts per electron at 20 keV, this relates to approximately 49.8 pA current on the detector, which is approximately 60 % of the incident current. A similar calculation on the Au/MoS₂ data indicates that 65% of the incident current was collected on the pixelated detector. The remainder would either be directed towards the HAADF detector (located above the pixelated detector and having a larger angular range), backscattered from the sample, or did not contribute enough energy to an individual pixel to pass the counting threshold of the detector, which is especially likely in the case of inelastic scattering with the specimen.

During the setup and alignment, the projector alignment coils were adjusted to give a distortion that visually appeared to be radial, which appeared coincident with the bright-field disc remaining approximately stationary while the projector lens was varied. Additionally, we aimed to position the bright-field disc away from the central sensor segment boundary which forms a cross-like artifact on the Quadro detector. As ptychographic algorithms require the unscattered beam to be at the centre of the diffraction images, the data was digitally shifted post-acquisition to recenter the beam in the diffraction plane. Regions outside to the original diffraction data were set to a value of zero. During reconstruction, these mask regions, which are outside the original field of collection and so do not represent real collected data, are left to float. The algorithm does not attempt to match the intensity in this extrapolated or set-to-zero region, or use this region in the evaluation of the error metric.
 
A ptychographic reconstruction was first made on the as collected Au/aC diffraction data. A radial average of the FT of the reconstruction was then fitted – using a pincushion distortion model as described later – to a powder diffraction pattern for gold, to determine the distortion present in the collected diffraction data. The extracted distortion model was subsequently used to undistort the diffraction data, prior to performing another reconstruction. Following this another comparison was made between the FT of the reconstruction and model diffraction data. If the model and reconstruction FT radial peak positions differed by more than 1%, another reconstruction and diffraction-model fit was performed, followed by a further distortion correction to the diffraction data. With our data, two reconstruction and fit cycles produced agreement between the peak positions of better than 0.5% out to the {135} diffraction ring of gold. Key parameters for the Au/aC reconstruction algorithm are given next, followed by a description of the diffraction data model fitting.
 
The methods of data collection are also described in the associated open access publication, linked to at https://github.com/ArthurBlackburn/PtychoRunner.

2. Methods for processing the data:

The processing of the data to form the ptychographic reconstructions here was done using a custom computational workflow based on the open-source PtychoShelves_EM MATLAB framework. We developed a wrapper tool, PtychoRunner, to manage batch processing and parameter control, and integrate additional analysis and constraint functions. We defined reconstruction parameters in spreadsheet files, which allows reproducible reconstructions and facilitates automated or parallel execution on GPU-enabled systems using job schedulers like SLURM.

Here the ptychographic reconstructions sequentially used the single-slice ePIE algorithm, an iterative least-squares maximum likelihood (LSQ-ML) solver for ptychography [4] implemented within the PtychoShelves code[5], and a multi-slice adaptation of (MS-)LSQ-ML, which was implemented within a modified version of PtychoShelves[6]. The reconstructions for the Au/aC (Au/MoS₂) samples used 50 iterations of ePIE followed by 300 iterations of LSQ-ML, and 600 (1700) iterations of MS-LSQ-ML. The first 300 iterations of LS-LSQ-ML used 2 slices, and the remaining iterations employed 4 slices separated by 20 Å, combined with 3 incoherently summed probe modes. To gain successful reconstructions in the MoS₂ sample, avoiding crystalline sample related artefacts in the reconstructed electron probe, we found it necessary to constrain the probe’s FT to be within an angular range less than the innermost extent of the first order diffraction discs. Given our probe-illumination angle, this constraint was reasonable and physical, simply confining the beam to have its known circular aperture constrained form. A similar constraint on the FT of the probe was applied during the Au/aC dataset reconstruction. This was found to improve the phyiscality of the reconstructed probe, though does not alter FRC resolution measure limit. 
 
To determine the distortion in the FT of the Au/aC reconstruction, we first determined a radial average of the FT’s squared magnitude. A distortion mapped ideal radial diffraction profile for randomly oriented gold powder, calculated using CrystalMaker (from CrystalMaker Software Limited, Oxfordshire, U.K.) with a 20 keV beam energy, was then fitted to this average profile by adjusting the overall scale or magnification, M, of this ideal diffraction intensity profile and introducing a pincushion distortion coefficient, k, defined in the relationship

	s = M*r*(1 + k*r^2) (1)
	
where s and r are the distorted and undistorted radial distances (in reciprocal space) of the diffraction intensity profile respectively. Here the fitting was performed by minimizing the Euclidean distance between the experimental and model radial profiles, using the simulated annealing algorithm, though we also found that genetic algorithms and pattern search methods yielded very similar results. 

The pincushion distortion, was corrected in our diffraction patterns using the distortion parameter, k, determined from an earlier Au/aC reconstruction. The magnification correction, given through M, was simply corrected for by adjusting the reconstructed pixel size input to the reconstruction algorithm, rather scaling the image. In the distortion correction we set M = 1. After distortion correction the reconstructed pixel size for the Au/aC and Au/MoS₂ dataset were 23.8 pm and 33.14 pm respectively. In the final 1000 iterations of the Au/MoS₂ reconstruction the diffraction data was padded by a factor of 2, resulting in a reconstructed pixel size of 16.57 pm using a super-resolution type approach. 
  
The most widely used image distortion correction routines effectively assume that pixels values represent delta function sampling of the image intensity, and consider where this point value must be mapped to. However, here the values of our diffraction data pixels represent the integrated diffracted wave intensity over a finite sized (75 μm square) square pixels. Thus, in the same manner as distortion correction has been applied to fMRI imaging[7] we adjusted the intensity of the pixel values in our distortion correction algorithm to take account of the compression (or, in general, expansion) of image patches between the distorted and undistorted spaces. This adjustment is generally not considered when correcting for general purpose image distortions, including for crystallographic purposes. However, for ptychography intensity adjustment we recognize this correction as important: without it, electron counts on the undistorted pattern would be modified in a non-linear manner, potentially corrupting the centre of mass of our diffraction discs. 

For example, imagine collecting pincushion distorted diffraction data where every pixel recorded one electron. Correcting the data without intensity adjustment would produce a smaller region of pixels with a value of one, surrounded by a region with an undefined value (or an arbitrarily assigned value, such as zero). Naively summing the pixel values in the ‘corrected’ image would then incorrectly indicate a smaller number of electron counts. Accordingly, prior to performing a conventional (non-intensity adjusting) distortion correction mapping we multiply the intensity of source data by

	 dA_s/dA_r = s*dtheta*ds/(r*dtheta*dr) = s*ds/(r*dr) = M^2*(1+k*r^2)*(1 + 3*k*r^2)	(2)
	 
In the manner consistent with that applied to fMRI image correction[7], the correction is determined by considering the ratio of infinitesimal areas of distorted and undistorted image space patches, dAs and dAr respectively. 
 
Prior works shows that a ptychographic reconstruction is only valid for a coherently scattering object thickness of T, given by the theoretical approximation,

	  T <= k*(δr)^2/λ,	(3)
	  
where δr is the image resolution, λ is the electron wavelength, and k is a constant in the range 2 - 5.2 [8]. With δr = 0.07 nm (our approximate resolution), and a beam energy of 20 keV giving λ = 8.59 pm, we thus should have T ≤ 1.1 – 3.0 nm for a single slice reconstruction approximation to hold. When the sample has a thickness greater than T, multi-slice ptychographic methods must be used (perhaps in combination tomography). Here, we opted to use n = 4 slices (where n is the number of slices), with their separation (T_s) approximately in the middle of range of T, thus setting T_s = 2 nm, to give a total acceptable object thickness of n*T_s = 8 nm.
 
Matching to our experimental diffraction patterns indicates that the MoS₂ specimen is composed of 5 layers (3.1 nm), which fits within the thickness limit given by n*T_s. For the gold islands, we estimate their thickness by taking their mean phase shift relative to the MoS₂ support from the experimental reconstructions and compare this to the phase shift seen in a reconstruction produced from simulated diffraction patterns for which the thickness of gold is known. Thus, we estimate our experimental gold island to have a thickness of 3.8 nm. This gives a total thickness of 6.9 nm for this region of the Au/MoS₂ reconstruction. Looking at an example Au particle from the Au/aC reconstruction, the peak phase shift difference between the particle and the amorphous carbon support is approximately 3 radians. A 20 keV electron beam has an interaction constant of 1.86 × 10^7 rad V^-1 m^-1, and gold has an estimated mean inner potential of 21.4 eV, thus implying a peak gold particle thickness of 7.5 nm, which is similar to the particle’s lateral dimensions, as would be expected for a near spherical particles.
 
Any hydrocarbon-based contamination on the sample is likely to be amorphous and mobile in nature and thus effectively contribute an incoherent background to the otherwise coherent diffraction data. Incoherent contributions to diffraction data cannot be effectively reconstructed by the ptychographic algorithms. Thus, the thickness limit consideration encapsulated in Equation (3) relates only to contiguous regions of the sample or object that scatter coherently: namely, the crystalline or poly-crystalline gold and MoS₂ in our samples.
 
Consequently, even though our samples likely had an amorphous carbon-based surface contamination layer (and in the case of the Au/aC sample a ~15 nm thickness amorphous carbon support film) that make the total sample thickness values greater than the 6.9 nm and 7.5 nm just determined, the reconstructions retain validity as the coherently scattering thickness regions are less than n*T_s (= 8 nm). However, the predominantly incoherent component to our diffraction data produced by the amorphous support or contamination material, likely contributes to intensity observed in the higher order probe model modes. This is more so in the absence of incoherent background estimation and removal, which was applied to the MoS₂ reconstruction but not to the Au/aC reconstruction.
 
Our simulated images and diffraction data for the Au – MoS₂ sample model were generated using Prismatic 2.0 software, using a multi-slice approach, [9, 10]. The single slice synthetic test object was prepared using Adobe Illustrator software combined with Matlab to overlay a series of double-gaussian peaks with variable separations at well separated locations throughout the image. The object was computationally ‘illuminated’ with a beam of the same energy (wavelength) as in our experiments, and with convergence half-angles, α, defoci, and overlap that were approximately in the same range as our experimental conditions. Specifically, α = 7.5, 10, 15 and 20 mrad were used in the simulations, where the latter two of these conditions would have resulted in overlap between the Bragg diffracted beams from their respective samples. These data were then reconstructed using the same parameters as for our experimental reconstruction of the Au-aC or Au-MoS₂ sample. The resulting reconstructed phase was then Fourier transformed and compared to the phase of our source object to determine a radially averaged signal to noise ratio (SNR), using the method in Eqn 44 of Thibault et al. [11], and phase contrast transfer function (PCTF).

The FRC was determined by splitting the ptychographic dataset into two subsets, forming the subsets from the even and odd numbered diffraction patterns (DP) as acquired, forming a chequerboard type pattern, respectively. Reconstructions were produced for each subset, using the same recipe as the full dataset reconstruction. These two reconstructions were effectively formed at half the electron dose of the full set. This method of producing independent subsets reduces the influence of sample drift and electromagnetic-noise induced beam-position movement upon on the FRC, as discussed and used in prior work[12, 13]. 

Effects and artefacts resulting from noise-induced random beam-position shifts upon individual reconstructions can be minimized or eliminated using correction algorithms[14, 15]. However, it is preferable to minimize or eliminate physical noise sources rather than computationally correct for them to reduce or remove the computational cost and residual uncertainty inherent to the correction algorithms. Thus, here we placed our SEM in a low-noise environment and observed stable, low-noise, high-resolution focussed probe conventional SE imaging prior to performing our defocused-probe data acquisitions. Applying positional correction algorithms within our reconstruction process resulted in insignificant changes of our Fourier error metric and did not improve our resolution measure. Thus, we did not apply such correction algorithms in the presented and analyzed reconstructions, and took the insignificant change of our Fourier error metric as evidence of low scan noise and good beam stability in our SEM. 

Each reconstruction produced summary images and high-resolution MATLAB data files that documented convergence behavior and reconstruction parameters.

The processing of the data is also described in the associated publication, and within the associate code provided at https://github.com/ArthurBlackburn/PtychoRunner. 

3. Instrument- or software-specific information needed to interpret the data:

The *.mat files above store data in Matlab v7.3 format. Version 7.3 MAT-files use an HDF5-based format that stores data in compressed chunks. The *.hp and *.hdf5 files above are HDF5 format, and can be opened using standard HDF5 libraries, such as described at https://github.com/HDFGroup/.

4. Standards and calibration information, if appropriate: 

The distortion corrected gold on amorhphous diffraction datasets have been calirbrated and corrected against the known lattice information for gold.


References

[1] Reidy, K. et al. Direct imaging and electronic structure modulation of moiré superlattices at the 2D/3D interface. Nature Communications 12, 1290 (2021). https://doi.org/10.1038/s41467-021-21363-5 
[2] Blackburn, A. M., Cordoba, C., Fitzpatrick, M. R. C. and Mcleod, R. A., 2025, Sub-ångström resolution ptychography in a scanning electron microscope at 20 keV, Submitted, Under Review. (see also https://github.com/ArthurBlackburn/PtychoRunner)
[3] Tinti, G. et al. The EIGER detector for low-energy electron microscopy and photoemission electron microscopy. Journal of Synchrotron Radiation 24, 963-974 (2017). https://doi.org/10.1107/S1600577517009109). 
[4] Odstrčil, M., Menzel, A. & Guizar-Sicairos, M. Iterative least-squares solver for generalized maximum-likelihood ptychography. Optics Express 26, 3108-3123 (2018). https://doi.org/10.1364/OE.26.003108
[5] Wakonig, K., Stadler, H.-C., Odstrcil, M., Tsai, E. H. R., Diaz, A., Holler, M., Usov, I., Raabe, J., Menzel, A. and Guizar-Sicairos, M., 2020, PtychoShelves, a versatile high-level framework for high-performance analysis of ptychographic data, Journal of Applied Crystallography 53, 574-586.
[6] Chen, Z. J., Yi; Muller, David A.; Odstrčil, Michal (2021) PtychoShelves_EM, source code for multislice electron ptychography, (Zenodo https://doi.org/10.5281/zenodo.4659690).
[7] Jezzard, P. & Balaban, R. S. Correction for geometric distortion in echo planar images from B0 field variations. Magnetic Resonance in Medicine 34, 65-73 (1995). https://doi.org/10.1002/mrm.1910340111
[8] Tsai, E. H. R., Usov, I., Diaz, A., Menzel, A. & Guizar-Sicairos, M. X-ray ptychography with extended depth of field. Optics Express 24, 29089-29108 (2016). https://doi.org/10.1364/OE.24.029089
[9] Rangel DaCosta, L. et al. Prismatic 2.0 – Simulation software for scanning and high resolution transmission electron microscopy (STEM and HRTEM). Micron 151, 103141 (2021). https://doi.org/10.1016/j.micron.2021.103141
[10] Fitzpatrick, M. R. C. Prismatique (2023).https://gitlab.com/mrfitzpa/prismatique
[11] Thibault, P. & Guizar-Sicairos, M. Maximum-likelihood refinement for coherent diffractive imaging. New Journal of Physics 14, 063004 (2012). https://doi.org/10.1088/1367-2630/14/6/063004
[12] Ding, Z. et al. Three-dimensional electron ptychography of organic–inorganic hybrid nanostructures. Nature Communications 13, 4787 (2022). https://doi.org/10.1038/s41467-022-32548-x
[13] Chen, Z. et al. Electron ptychography achieves atomic-resolution limits set by lattice vibrations. Science 372, 826 (2021). https://doi.org/10.1126/science.abg2533
[14] Dwivedi, P., Konijnenberg, A. P., Pereira, S. F. & Urbach, H. P. Lateral position correction in ptychography using the gradient of intensity patterns. Ultramicroscopy 192, 29-36 (2018). https://doi.org/10.1016/j.ultramic.2018.04.004
[15] Maiden, A. M., Humphry, M. J., Sarahan, M. C., Kraus, B. & Rodenburg, J. M. An annealing algorithm to correct positioning errors in ptychography. Ultramicroscopy 120, 64-72 (2012). https://doi.org/10.1016/j.ultramic.2012.06.001

