function rgbCmap = GiveMeColors(whatColors)

switch whatColors
case 'areas'
    rgbCmap = [57,146,131;
                212,96,82;
                102,220,169;
                188,148,149;
                161,207,207;
                196,206,80]/255;
otherwise
    error('Unknown map: ''%s''',whatColors);
end
