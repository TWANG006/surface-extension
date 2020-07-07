function [Xa, Ya, Za] = CleanNaNFromSurfData(X,Y,Z)
% Function
%   CleanCircSurfData(X,Y,Z)
% Purpose
%   Clean the Circular surface map to include only the valid data points
% Inputs:
%   X, Y, Z: coordinates & surf data before cleaning
% Outputs:
%   Xa, Ya, Za: coordinates & surf data After cleaning   

ids_nanZ = ~isnan(Z);

pys = find(Y(:,1)==nanmin(nanmin(Y(ids_nanZ))));
pye = find(Y(:,1)==nanmax(nanmax(Y(ids_nanZ))));

pxs = find(X(1,:)==nanmin(nanmin(X(ids_nanZ))));
pxe = find(X(1,:)==nanmax(nanmax(X(ids_nanZ))));

Xa = X(pys:pye, pxs:pxe);
Ya = Y(pys:pye, pxs:pxe);
Za = Z(pys:pye, pxs:pxe);

end