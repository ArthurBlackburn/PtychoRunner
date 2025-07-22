function a1String = idx2A1(idx)
% This function translates an index to the Excel A1 notation that is used in xlswrite.
% 
% Matt Brunner (2023). Convert index to Excel A1 notation 
% (https://www.mathworks.com/matlabcentral/fileexchange/
% 28794-convert-index-to-excel-a1-notation), 
% MATLAB Central File Exchange

alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

if idx < 27
    a1String = alphabet(idx);
else
    idx2 = rem(idx,26);
    if idx2 == 0
        a1String = [alphabet(floor(idx/26)-1),'Z'];
    else
        a1String = [alphabet(floor(idx/26)),alphabet(idx2)];
    end
end

end
