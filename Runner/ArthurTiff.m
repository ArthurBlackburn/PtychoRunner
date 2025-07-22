function ArthurTiff(data_in, filename, description_ascii, unit_desc)
% Just a simple function for saving 2-d and 3-d matrices as a TIF file. As TIF supports 64-bit 
% (double) pixel values, uses lossless compression, and is read by many programs it is a very useful 
% format. That said, not all programs can read in 64 bit values. Photoshop etc, does a good job - 
% Windows browser does not. The bit depth is set by the data type of the input. Look inside the 
% function!
%
% Usage:
% ArthurTiff(data_in, filename, description_ascii, unit_desc).
%
% data_in:  The data in!
% filename: The filename to save to. If no filename supplied the function does nothing, but issues a
%           warning.
% description_ascii: Tif files can have an ascii based description addded. Add this here as string
%                    or char array, or for example serialized JSON info.
% unit_desc: Specify the pixel size and unit name in a cell array. E.g. if the pixels are
%            20 nm, then use {20, 'nm'}. Internally the file uses Tiff.ResolutionUnit.Centimeter and
%            conversions are made. See the code internals...
%
% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2018
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
%
% *******************************************************************



if isempty(filename)
    warning('No filename supplied, not saving image.');
    return
end

if nargin < 3
     description_ascii = [];
end

if nargin < 4
    unit_desc = [];
end

% Construct TIFF image...
t = Tiff(filename, 'w'); 

% ...with these default parameters...
tagstruct = struct(...
    'ImageLength'        , size(data_in,1),...
    'ImageWidth'         , size(data_in,2),...
    'Compression'        , Tiff.Compression.LZW,...%     
    'SampleFormat'       , Tiff.SampleFormat.IEEEFP,...  % floating point
    'Photometric'        , Tiff.Photometric.MinIsBlack,... %   
    'BitsPerSample'      , 32,... % 8 bytes / double
    'SamplesPerPixel'    , 1,...
    'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);

% ... and with the bit depth changed accordinf to the dtype ...
dtype = class(data_in);
switch dtype
    case 'double'
        tagstruct.SampleFormat =  Tiff.SampleFormat.IEEEFP;
        tagstruct.BitsPerSample = 64; 
        % though beware that not many programs can read 64-bit tiff, and
        % that in mnay cases it is overkill for display purposes.        
    case 'single'
        tagstruct.SampleFormat =  Tiff.SampleFormat.IEEEFP;
        tagstruct.BitsPerSample = 32;         
    case 'uint32'
        tagstruct.SampleFormat =  Tiff.SampleFormat.UInt;
        tagstruct.BitsPerSample = 32;         
    case 'uint16'
        tagstruct.SampleFormat =  Tiff.SampleFormat.UInt;
        tagstruct.BitsPerSample = 16;
    case 'uint8'
        tagstruct.SampleFormat =  Tiff.SampleFormat.UInt;
        tagstruct.BitsPerSample = 8;         
end

if ~isempty(description_ascii)
    tagstruct.ImageDescription = description_ascii;
end
if iscell(unit_desc)
    tagstruct.ResolutionUnit = Tiff.ResolutionUnit.Centimeter;
    if ischar(unit_desc{2})
        im_unit = unit_desc{2};
        switch im_unit
            case 'm'
                multip = 100;
            case 'inch'
                multip = 2.54;
            case 'cm'
                multip = 1;
            case 'mm'
                multip = 1e-1;
            case 'nm'
                multip = 1e-7;
            case 'µm'
                multip = 1e-4;
            case 'um'
                multip = 1e-4;
            case 'micron'
                multip = 1e-4;
            case 'pm'
                multip = 1e-10;
            case 'A'
                multip = 1e-8;
            otherwise
                multip = 1;
                warning('Unknown unit specified, using default centimeters');
        end
    end
    tagstruct.XResolution = unit_desc{1}*multip;
    tagstruct.YResolution = unit_desc{1}*multip;
end

depth = size(data_in,3);

for d = 1:depth
    t.setTag(tagstruct);
    t.write(data_in(:, :, d));
    if d ~= depth
       t.writeDirectory();
    end
end

t.close();