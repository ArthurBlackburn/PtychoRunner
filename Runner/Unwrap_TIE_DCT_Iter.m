%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D phase Unwrapping algorithm based on a manuscript entitled "Robust 2D phase unwrapping algorithm based
% on the transport of intensity equation",which was submitted to Measurement Science and Technology(MST).
% Inputs:
%   * phase_wrap: wrapped phase from -pi to pi
% Output:
%   * phase_unwrap: unwrapped phase 
%   * N: number of iterations 
% Author:Zixin Zhao (Xi'an Jiaotong University, 08-15-2018)
% Email:zixinzhao@xjtu.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phase_unwrap,N]=Unwrap_TIE_DCT_Iter(phase_wrap, timeout, maxIters)
   % make sure we have the right version of dct.
   % this is where we get laid bare just how bad Matlab can be...
   % As this is the official workaround...
   olddir = cd(fullfile(toolboxdir('signal'),'signal')); % use this one!
%    olddir = cd(fullfile(toolboxdir('images'),'images','private'));
%    REQD_dct = @(x, varargin) dct(x, varargin{:});
%    REQD_idct = @(b, varargin) idct(b,varargin{:});
   REQD_dct = str2func('dct');
   REQD_idct = str2func('idct');

   cd(olddir);

   if nargin < 2
       timeout = Inf;
   end
   
   if nargin < 3
       maxIters = Inf;
   end
   
   phi1 = unwrap_TIE(phase_wrap, REQD_dct, REQD_idct );
   phi1=phi1+mean2(phase_wrap)-mean2(phi1); %adjust piston
    K1=round((phi1-phase_wrap)/2/pi);  %calculate integer K
    phase_unwrap=phase_wrap+2*K1*pi; 
    residue=wrapToPi(phase_unwrap-phi1);
    phi1=phi1+unwrap_TIE(residue, REQD_dct, REQD_idct );
    phi1=phi1+mean2(phase_wrap)-mean2(phi1); %adjust piston
    K2=round((phi1-phase_wrap)/2/pi);  %calculate integer K
    phase_unwrap=phase_wrap+2*K2*pi; 
    residue=wrapToPi(phase_unwrap-phi1);
    N=0;
   tic;
   while ((sum(sum(abs(K2-K1)))>0) && (toc < timeout) && N < maxIters)
       K1=K2;
       phic=unwrap_TIE(residue, REQD_dct, REQD_idct );
     phi1=phi1+phic;
     phi1=phi1+mean2(phase_wrap)-mean2(phi1); %adjust piston
    K2=round((phi1-phase_wrap)/2/pi);  %calculate integer K
    phase_unwrap=phase_wrap+2*K2*pi; 
    residue=wrapToPi(phase_unwrap-phi1);
    N=N+1;
   end
   if toc > timeout
       fprintf('NOTE: Unwrap_TIE_DCT_Iter terminated by timeout of %1.0f sec\n',timeout);
   end
end

function [phase_unwrap]=unwrap_TIE(phase_wrap, dct, idct)
      psi=exp(1i*phase_wrap);
%       edx = [zeros([size(psi,1),1]), wrapToPi(diff(psi, 1, 2)), zeros([size(psi,1),1])]; % original code
%       edy = [zeros([1,size(psi,2)]); wrapToPi(diff(psi, 1, 1)); zeros([1,size(psi,2)])]; % original code
     % wrapToPi on a complex number (diff of psi), is 'unusual', as wraptoPi operates on real numbers.
     % even if say wrap to pi was modified to operate on complex number
     % arg, that would be a pointless operation as operation as diff of 2 complex
     % numbers are always going to have phase -pi<= phase < pi...
      edx = [zeros([size(psi,1),1]), diff(psi, 1, 2), zeros([size(psi,1),1])];
      edy = [zeros([1,size(psi,2)]); diff(psi, 1, 1); zeros([1,size(psi,2)])];
      lap = diff(edx, 1, 2) + diff(edy, 1, 1); %calculate Laplacian using the finite difference
      rho=imag(conj(psi).*lap);   % calculate right hand side of Eq.(4) in the manuscript
   phase_unwrap = solvePoisson(rho, dct, idct); 
end

function phi = solvePoisson(rho, dct, idct)
    % solve the poisson equation using DCT
    dctRho = dct2(rho, [], [], dct);
    [N, M] = size(rho);
    [I, J] = meshgrid([0:M-1], [0:N-1]);
    dctPhi = dctRho ./ 2 ./ (cos(pi*I/M) + cos(pi*J/N) - 2);
    dctPhi(1,1) = 0; % handling the inf/nan value
    % now invert to get the result
    phi = idct2(dctPhi, [], [], idct);
end

% Put the below here, otherwise dct and icdt in the image toolbox is used, which
% then calls dct / idct in the image toolbox rather than the one in the signal
% toolbox which can handle gpuArrays. Be sure to put signal toolbox above
% image toolbox, or use the above 'hack above' to pass handle to correct
% dct function in signal toolbox that can handle gpuArrays.

function b=dct2(arg1,mrows,ncols, dct)

[m, n] = size(arg1);
% Basic algorithm.
if (nargin == 1) || isempty(mrows)
  if (m > 1) && (n > 1)
    b = dct(dct(arg1).').';
    return;
  else
    mrows = m;
    ncols = n;
  end
end

% Padding for vector input.
a = arg1;
if nargin==2 || isempty(ncols)
    ncols = mrows(2);
    mrows = mrows(1); 
end
mpad = mrows; npad = ncols;
if m == 1 && mpad > m, a(2, 1) = 0; m = 2; end
if n == 1 && npad > n, a(1, 2) = 0; n = 2; end
if m == 1, mpad = npad; npad = 1; end   % For row vector.

% Transform.

b = dct(a, mpad);
if m > 1 && n > 1, b = dct(b.', npad).'; end

end

function a = idct2(arg1,mrows,ncols, idct)
%IDCT2 2-D inverse discrete cosine transform.
%   B = IDCT2(A) returns the two-dimensional inverse discrete
%   cosine transform of A.
%
%   B = IDCT2(A,[M N]) or B = IDCT2(A,M,N) pads A with zeros (or
%   truncates A) to create a matrix of size M-by-N before
%   transforming. 
%
%   For any A, IDCT2(DCT2(A)) equals A to within roundoff error.
%
%   The discrete cosine transform is often used for image
%   compression applications.
%
%   Class Support
%   -------------
%   The input matrix A can be of class double or of any
%   numeric class. The output matrix B is of class double.
%
%   References:
%   -----------
%   [1] A. K. Jain, "Fundamentals of Digital Image
%       Processing", pp. 150-153.
%   [2] Wallace, "The JPEG Still Picture Compression Standard",
%       Communications of the ACM, April 1991.
%
%   Example
%   -------
%   % Perform inverse DCT on an image.
%       RGB = imread('autumn.tif');
%       I = rgb2gray(RGB);
%       J = dct2(I);
%       imshow(log(abs(J)),[]), colormap(gca,jet), colorbar
%
%   % The commands below set values less than magnitude 10 in the
%   % DCT matrix to zero, then reconstruct the image using the
%   % inverse DCT function IDCT2.
%
%       J(abs(J)<10) = 0;
%       K = idct2(J);
%       figure, imshow(I)
%       figure, imshow(K,[0 255])
%
%   See also DCT2, DCTMTX, FFT2, IFFT2.

%   Copyright 1992-2018 The MathWorks, Inc.



[m, n] = size(arg1);
% Basic algorithm.
if (nargin == 1) || isempty(mrows)
  if (m > 1) && (n > 1)
    a = idct(idct(arg1).').';
    return;
  else
    mrows = m;
    ncols = n;
  end
end

% Padding for vector input.

b = arg1;
if nargin==2 || isempty(ncols)
    ncols = mrows(2); 
    mrows = mrows(1); 
end

mpad = mrows; npad = ncols;
if m == 1 && mpad > m, b(2, 1) = 0; m = 2; end
if n == 1 && npad > n, b(1, 2) = 0; n = 2; end
if m == 1, mpad = npad; npad = 1; end   % For row vector.

% Transform.

a = idct(b, mpad);
if m > 1 && n > 1, a = idct(a.', npad).'; end

end

