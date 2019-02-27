function LabelCurrentAxes(theText,ax,fontSize,whereToAnnotate,fontColor)
% ------------------------------------------------------------------------------
% Labels a subplot with 'a', 'b', etc.
% Ben Fulcher, 2015-02-26
% ------------------------------------------------------------------------------

if nargin < 2
    ax = gca;
end
if nargin < 3
    fontSize = 18;
end
if nargin < 4
    whereToAnnotate = 'topRight';
end
if nargin < 5
    fontColor = 'k';
end

% ------------------------------------------------------------------------------
switch whereToAnnotate
case 'topRight'
    x = ax.XLim(1) + (ax.XLim(2)-ax.XLim(1))*0.9;
    y = ax.YLim(1) + (ax.YLim(2)-ax.YLim(1))*0.9;
case 'topLeft'
    x = ax.XLim(1) + (ax.XLim(2)-ax.XLim(1))*0.1;
    y = ax.YLim(1) + (ax.YLim(2)-ax.YLim(1))*0.9;
otherwise
    error('Unknown annotation position');
end

text(ax,x,y,theText,'FontSize',fontSize,'FontWeight','bold','Color',fontColor);

end
