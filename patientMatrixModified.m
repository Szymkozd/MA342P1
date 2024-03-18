function [k] = patientMatrixModified(x,y,s,oar,s_oar,plt)
%PATIENTMATRIX computes the patient matrix
%inputs:     x: x coordinate of center of tumor square
%            y: y coordinate of center of tumor square
%            s: size length of tumor
%          oar: organs at risk matrix (2xn, each col is coordinates of an OAR)
%        s_oar: size of each organ at risk
%          plt: whether to plot the resulting matrix (1) or not (0)
%
%outputs:    k: matrix encoding the following
%                  0: outside patient (center of square is outside)
%                  1: normal tissue
%                  2: tumor
% 
% grid is given by -50 <= x,y <= 50

k = zeros(500,500); % initialize all points to be outside matrix

% mark points inside circle
for a = 201:300 % row
    for b = 201:300 % col
        realC = matIdxToCoords(a,b);
        if (norm(realC) < 50) % patient radius is 10 cm = 50 pixels
            k(a,b) = 1;
        end
    end
end

% mark tumor
[r,c] = coordsToMatIdx(x,y);
for a = round(r-s/2):round(r+s/2) % row
    for b = round(c-s/2):round(c+s/2) % col
        k(a,b) = 2;
    end
end

% mark OARs
for i = 1:size(oar,2)
    [r,c] = coordsToMatIdx(oar(1, i), oar(2,i));
    for a = round(r-s_oar/2):round(r+s_oar/2) % row
        for b = round(c-s_oar/2):round(c+s_oar/2) % col
            k(a,b) = 3;
        end
    end
end

if (plt == 1)
    heatmap(k,'GridVisible','off');
end

    function [v] = matIdxToCoords(r,c) 
        v = [c-250.5, 250.5-r];
    end

    function [r,c] = coordsToMatIdx(x,y)
        r = 250.5 - y;
        c = x + 250.5;
    end
end