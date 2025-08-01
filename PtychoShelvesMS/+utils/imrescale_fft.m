% IMRESCALE_FFT subpixel precision rescaling based on multiplication by a 
% matrix of fourier transformation, fast only for small arrays 
% Inputs:  
%     **img   - 2D or stack of 2D images 
%     **scale - scaling factor 
% *returns*:
%     ++img   - 2D or stack of 2D images scaled by factor scale

%*-----------------------------------------------------------------------*
%|                                                                       |
%|  Except where otherwise noted, this work is licensed under a          |
%|  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
%|  International (CC BY-NC-SA 4.0) license.                             |
%|                                                                       |
%|  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
%|                                                                       |
%|       Author: CXS group, PSI                                          |
%*-----------------------------------------------------------------------*
% You may use this code with the following provisions:
%
% If the code is fully or partially redistributed, or rewritten in another
%   computing language this notice should be included in the redistribution.
%
% If this code, or subfunctions or parts of it, is used for research in a 
%   publication or if it is fully or partially rewritten for another 
%   computing language the authors and institution should be acknowledged 
%   in written form in the publication: “Data processing was carried out 
%   using the “cSAXS matlab package” developed by the CXS group,
%   Paul Scherrer Institut, Switzerland.” 
%   Variations on the latter text can be incorporated upon discussion with 
%   the CXS group if needed to more specifically reflect the use of the package 
%   for the published work.
%
% A publication that focuses on describing features, or parameters, that
%    are already existing in the code should be first discussed with the
%    authors.
%   
% This code and subroutines are part of a continuous development, they 
%    are provided “as they are” without guarantees or liability on part
%    of PSI or the authors. It is the user responsibility to ensure its 
%    proper use and the correctness of the results.



function img_rescale = imrescale_fft(img, scale)

    if scale == 1 || isnan(scale)
        img_rescale = img;
        return
    end

    
    N = size(img);

    for i = 1:2
        if N(1) ~= N(2) || i == 1
            ind_x = real(zeros(N(i),1, 'like', img));
            ind_x(:) = (0:N(i)-1)/N(i)-0.5; 
            grid = ind_x*ind_x';
            grid = -2i*pi*N(i)/scale*grid;
            W{i} = exp(grid)'/N(i);  % matrix of fourier transformation 
        else
            W{2} = W{1};
        end
    end    
    
    fimg = fftshift(fft2(fftshift(img)));

    
    if size(img,3) > 1
        fimg2 = W{1}*reshape(fimg,N(1),[]);
        fimg2 = reshape(fimg2, N(1),N(2),[]);
        fimg2 = permute(fimg2, [2,1,3]);
        fimg2 = reshape(fimg2, N(2),[]);
        img_rescale = (W{1}*fimg2);  % rescale and fft back 
        img_rescale = reshape(img_rescale, N(2),N(1),[]);
        img_rescale = permute(img_rescale, [2,1,3]);
    else
        fimg2 = W{1}*fimg;
        img_rescale = (W{2}*fimg2.').';
    end
    
end

