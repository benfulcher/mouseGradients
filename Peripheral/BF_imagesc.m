function BF_imagesc(X,whatColorNaN)
% Plots an imagesc, and then adds black rectangles over NaN values

if nargin < 2
    whatColorNaN = 'k';
end

%-------------------------------------------------------------------------------
% Plot the color matrix:
%-------------------------------------------------------------------------------
imagesc(X);

%-------------------------------------------------------------------------------
% Add rectangles over NaN values:
%-------------------------------------------------------------------------------

if any(isnan(X(:)))
    [theNaNs_i,theNaNs_j] = find(isnan(X));
    for i = 1:length(theNaNs_i)
        rectangle('Position',[theNaNs_j(i)-0.5,theNaNs_i(i)-0.5,1,1],...
                        'FaceColor',whatColorNaN, ...
                        'EdgeColor',whatColorNaN)
    end
end


end
