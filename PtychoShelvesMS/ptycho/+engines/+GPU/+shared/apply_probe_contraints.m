% APPLY_PROBE_CONTRAINTS apply support constrains on the probe in the real space, fourier space or any other
% plane if provided, it useds factor in mode.support_back_propagation_factor to perform ASM propagation 
%
% probe = apply_probe_contraints(probe, mode)
%
% ** probe     complex array with probe / probes 
% ** mode       structure containing parameters for selected probe mode 
%
% returns:
% ** probe     complex array with probe / probes 
%
%

% Academic License Agreement
%
% Source Code
%
% Introduction 
% �	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE�s responsibility to ensure its proper use and the correctness of the results.�
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379�382 (2008). 
%   (doi: 10.1126/science.1158573),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68�71 (2013). (doi: 10.1038/nature11806),
% for LSQ-ML method 
% M. Odstrcil, A. Menzel, M.G. Sicairos,  Iterative least-squares solver for generalized maximum-likelihood ptychography, Optics Express, 2018
% for OPRP method 
%  M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  "Ptychographic coherent diffractive imaging with orthogonal probe relaxation." Optics express 24.8 (2016): 8360-8369
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089�29108 (2016). 
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       � All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Z�rich, Switzerland. 


function probe = apply_probe_contraints(probe, mode)
    import math.*
    import utils.*
    import engines.GPU.shared.*

    if  ~isempty(mode.probe_support)
        % apply support contraint in real space (if nearfield propagated )
        if ~isempty(mode.support_fwd_propagation_factor)
            if isscalar(mode.support_fwd_propagation_factor) && isinf(mode.support_fwd_propagation_factor)
                probe =  fftshift_2D(fft2(fftshift_2D(probe)));  % propagate to infinity
            else
                probe =  ifft2(fft2(probe) .* mode.support_propagation_factor);     
            end
        end   

        %% apply real-space support 
        probe = probe .* mode.probe_support;

        if ~isempty(mode.support_back_propagation_factor)
            if isscalar(mode.support_back_propagation_factor) && isinf(mode.support_back_propagation_factor)
                probe =  fftshift_2D(ifft2(fftshift_2D(probe)));  % propagate to infinity
            else
                probe =  ifft2(fft2(probe).* mode.support_back_propagation_factor);  
             end
        end   
    
    end
    
    Np_p = size(probe); 
    verbose(2,'For debug: mode.probe_scale_upd = %3.3e',mode.probe_scale_upd(end))
%     if mode.probe_scale_upd(end) > 0  && ~isempty(mode.probe_scale_window)
% above is original; below is modified
    if mode.probe_scale_upd(end) < 0  && ~isempty(mode.probe_scale_window)
        % apply windowing to avoid boundary issues when subpixel probe
        % rescaling is used.
        % Don't think is actually needed, as scaling sets all outside of a
        % reduced scale probe to zero anyway... but for now let's leave it
        % here. Also seems to me the inequaliy was initially the wrong way
        % around. -Ve scale produces smaller probe in old map area.
        verbose(2,'For debug: applied scale window on probe with mode.probe_scale_upd = %3.3e',mode.probe_scale_upd(end))
        probe = probe .* mode.probe_scale_window ; % here we are apply probe_scale window on probe in real space
        % at this point the probe is not fftshifted.
    end
    if  ~isempty(mode.probe_support_fft) || mode.probe_scale_upd(end) ~= 0 
        %% apply contraint in the detector plane
        
        % propagate probe on the detector 
        probe = fwd_fourier_proj(probe, mode);  
        % now the probe is projected to the detector, i.e. fft
                    
        if ~isempty(mode.probe_support_fft)
            probe = probe .* mode.probe_support_fft;
        end
%         if mode.probe_scale_upd(end) < 0 && ~isempty(mode.probe_scale_window)
%             probe = probe .*  mode.probe_scale_window ; % here the same window on the probe but now in fft space... why?
%             % at the least it shoud be fftshifted, as was partly done in
%             % ptychosolver (but previously not consistently, being
%             % fftshifted when <0 but not when >0). But at this moment I see
%             % no reason to be apply an apodizing mask on image and fft
%             % plane, and on when scale_upd is either >0 or <0 from
%             % conditional above this one. Arthur B. 
%         end
        
        % propagate probe back to the sample plane 
        probe = back_fourier_proj(probe, mode);  
    end    
end

