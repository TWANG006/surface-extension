function [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension_Polyfit(...
    X, Y, Z,...unextended surface error map
    tif_mpp,...TIF sampling interval [m/pxl]
    Z_tif,...TIF profile
    order_m, order_n,...polynomial orders in y, x
    type...Chebyshev or Legendre
    )
% Function
%   [X_ext, Y_ext, Z_ext] = Surface_Extension_Polyfit(X, Y, Z,...
%                           X_tif, Y_tif, Z_tif, order_m, order_n, type)
% Purpose
%   Extend the surface error map using polynomial fitting

%% 0. Obtain required parameters
% Sampling intervals
surf_mpp = median(diff(X(1,:)));    ... surface sampling interval [m/pxl]

m = size(Z,1);  ... CA height [pixel]
n = size(Z,2);  ... CA width [pixel]

m_ext = round(tif_mpp*(size(Z_tif, 1))*0.5/surf_mpp);   ... extension size in y [pixel]
n_ext = round(tif_mpp*(size(Z_tif, 2))*0.5/surf_mpp);   ... extension size in x [pixel] 

ca_range.v_s = m_ext + 1;   ca_range.v_e = ca_range.v_s + m - 1;   ... y start & end ids of CA in FA [pixel]
ca_range.u_s = n_ext + 1;   ca_range.u_e = ca_range.u_s + n - 1;   ... x start & end ids of CA in FA [pixel]


%% 1. Initial extension matrices
% extension sizes
[X_ext, Y_ext] = meshgrid(-n_ext:n-1+n_ext, -m_ext:m-1+m_ext);  ...extension grid
X_ext = X_ext * surf_mpp + X(1,1);  ... adjust X grid add X(1,1)
Y_ext = Y_ext * surf_mpp + Y(1,1);  ... adjust Y grid add Y(1,1)
Z_ext = NaN(size(X_ext));   ... mark the Z_ext to NaN
Z_ext(ca_range.v_s:ca_range.v_e, ca_range.u_s:ca_range.u_e) = Z;... fill in the valid data point

%% Fit the edge values
% zero boundary contidion
w = 100;
Z_ext(:, 1) = 0;
Z_ext(1, :) = 0;
Z_ext(:, end) = 0;
Z_ext(end, :) = 0;

W = ones(size(Z_ext));
W(:, 1) = w;
W(1, :) = w;
W(:, end) = w;
W(end, :) = w;

%% 2. poly fit
[p, q] = meshgrid(0:order_n, 0:order_m);
X_nor = -1 + 2.*(X_ext - min(X_ext(:)))./(max(X_ext(:)) - min(X_ext(:)));
Y_nor = -1 + 2.*(Y_ext - min(Y_ext(:)))./(max(Y_ext(:)) - min(Y_ext(:)));

if(strcmp(type,'Chebyshev'))
    [z3, ~, ~] = ChebyshevXYnm(X_nor, Y_nor, p(:), q(:));
elseif(strcmp(type,'Legendre'))
    [z3, ~, ~] = LegendreXYnm(X_nor, Y_nor, p(:), q(:));
else
    error('Unkown polynomial type.');
end

z3_res = reshape(z3, [],size(z3,3));

A = z3_res(~isnan(Z_ext(:)),:);
b = Z_ext(~isnan(Z_ext(:)));

c = lscov(A,b, W(~isnan(Z_ext(:))));

for i = 1:length(c)
    z3(:,:,i) = z3(:,:,i)*c(i);
end

Z_ext = sum(z3,3);
% Z_ext(ca_range.v_s:ca_range.v_e, ca_range.u_s:ca_range.u_e) = Z;

end