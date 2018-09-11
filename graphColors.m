function [colors, colorNames] = graphColors(nColors,graphType)

% Create colors for graphs based on Maureen Stone's work for best colors to
% use on figures

% Inputs:
% nColors =  how many different colors needed?
% graphType = 1 for bar, 0 for lines and points

% Outputs:
% colors = an nColors x 3 (rgb) matrix of colors scaled from 0 to 1
% colorNames = an nColors x 1 cell array with the associated color names

switch graphType
    case 1 % bar
        blue = [114 147 203];
        orange = [225 151 76];
        green = [132 186 91];
        red = [211 94 96];
        gray = [128 133 133];
        purple = [144 103 167];
        brown = [171 104 194];
        gold = [204 194 16];
    otherwise % lines/points
        blue = [157 106 177];
        orange = [218 124 48];
        green = [62 150 81];
        red = [204 37 41];
        gray = [83 81 84];
        purple = [107 76 154];
        brown = [146 36 40];
        gold = [148 139 61];
end

colorMatrix = [blue;orange;green;red;gray;purple;brown;gold];
colorNames = {'blue';'orange';'green';'red';'gray';'purple';'brown';'gold'};

if nColors > 8
    nColors = 8;
end


colors = colorMatrix(1:nColors,:)/255;