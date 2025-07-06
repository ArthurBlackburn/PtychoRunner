# Codes for realization of sub-ångström resolution ptychography in a scanning electron microscope at 20 keV

Arthur M. Blackburn, University of Victoria, 2025

## Introduction

The program codes given here produce the key results of our work related to the realization of sub-ångström resolution ptychography in a scanning electron microscope at 20 keV, which have recently been submitted for publication[1]. The code here provides a wrapper for and modifies a version of Ptychoshelves_EM  [2, 3], which in turn is a modified subset of a prior release of PtychoShelves[4]. 

To allow some traceability back to the Ptychoshelves_EM  [2, 3] code, the first commit to the repository here is a copy of the code hosted on Zenodo [3], with some superfluous backup (\*.bak,  \*.asv) files removed, and some – but not all - Matlab based Unix / Windows mixed line-endings resolved. Modifications we have made to the Ptychoshelves_EM code can be seen by viewing the repository changelog.

A set of wrapper codes and modifications, called 'Ptychorunner' will subsequently be provided herewith. These facilitate using PtychoShelves_EM and PtychoShelves in a batch run environment, where all the parameters controlling the reconstructions are contained within \*.xls spreadsheet files. 

## References

1. Blackburn, A. M., Cordoba, C., Fitzpatrick, M. R. C. and Mcleod, R. A., 2025,  Sub-ångström resolution ptychography in a scanning electron microscope at 20 keV, Nature Communications, Under Review.
2. Chen, Z., Jiang, Y., Shao, Y.-T., Holtz, M. E., Odstrčil, M., Guizar-Sicairos, M., Hanke, I., Ganschow, S., Schlom, D. G. and Muller, D. A., 2021,  Electron ptychography achieves atomic-resolution limits set by lattice vibrations, Science 372, 826.
3. Chen, Z. J., Yi; Muller, David A.; Odstrčil, Michal (2021) PtychoShelves_EM, source code for multislice electron ptychography, (Zenodo  https://doi.org/10.5281/zenodo.4659690 ).
4. Wakonig, K., Stadler, H.-C., Odstrcil, M., Tsai, E. H. R., Diaz, A., Holler, M., Usov, I., Raabe, J., Menzel, A. and Guizar-Sicairos, M., 2020,  PtychoShelves, a versatile high-level framework for high-performance analysis of ptychographic dataThis article will form part of a virtual special issue of the journal on ptychography software and technical developments, Journal of Applied Crystallography 53, 574-586.
