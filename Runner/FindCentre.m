function centrePoint = FindCentre(input, show_op)
%% FindCentre()
%  Find the centre of mass of an image
%
% Inputs
%  input - An image or array of which to find the centre of mass.
%
% Outputs
%  centrePoint - Co-ordinates of the centre of mass.

    if nargin < 2
        show_op = false;
    end
    input = abs(input);
    [rowSize,colSize] = size(input);
     
    oneD = rowSize == 1 | colSize == 1;
    
    if oneD
        centrePoint = length(input)*(1 - mean((cumsum(input)/sum(input)))) + 1;
    else
        
        colCentre = colSize*(1 - mean((cumsum(sum(input))/sum(sum(input))))) + 1;
        rowCentre = rowSize*(1 - mean((cumsum(sum(input.'))/sum(sum(input.'))))) + 1;
        
        centrePoint = [rowCentre colCentre]; % i.e. y, x
        if show_op
            imagesc(input);
            line(centrePoint(2),centrePoint(1),'marker','x','color','y');
        end
    end
end
    
    
