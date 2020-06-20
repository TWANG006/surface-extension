function [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension_GP(...
    X, Y, Z,...unextended surface error map
    tif_mpp,...TIF sampling interval [m/pxl]
    Z_tif,...TIF profile
    fx_range,...
    fy_range...
    )
%% 0. Obtain required parameters
% Sampling intervals
surf_mpp = median(diff(X(1,:)));    ... surface sampling interval [m/pxl]

m = size(Z,1);  ... CA height [pixel]
n = size(Z,2);  ... CA width [pixel]

m_ext = round(tif_mpp*(size(Z_tif, 1))*0.5/surf_mpp);   ... extension size in y [pixel]
n_ext = round(tif_mpp*(size(Z_tif, 2))*0.5/surf_mpp);   ... extension size in x [pixel] 

ca_range.y_s = m_ext + 1;   ca_range.y_e = ca_range.y_s + m - 1;   ... y start & end ids of CA in FA [pixel]
ca_range.x_s = n_ext + 1;   ca_range.x_e = ca_range.x_s + n - 1;   ... x start & end ids of CA in FA [pixel]


%% 1. Initial extension matrices
% extension sizes
[X_ext, Y_ext] = meshgrid(-n_ext:n-1+n_ext, -m_ext:m-1+m_ext);  ...extension grid
X_ext = X_ext * surf_mpp + X(1,1);  ... adjust X grid add X(1,1)
Y_ext = Y_ext * surf_mpp + Y(1,1);  ... adjust Y grid add Y(1,1)

Z_ini = NaN(size(X_ext));
Z_ini(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e) = Z;... fill in the valid data point

Z_ext = zeros(size(X_ext));   ... mark the Z_ext to NaN
Z_ext(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e) = Z;... fill in the valid data point
Z_ext(isnan(Z_ext)) = 0;

G = 0*Z_ext; Gy = G;
G(~isnan(Z_ini)) = 1;... fill in the valid data point
Gy(ca_range.y_s:ca_range.y_e,:) = 1;
Gox = 0*Z_ext;
Goy = 0*Z_ext;
Gox(:, fix(size(Z_ext,2)/2)+1+fx_range)=1;
Goy(fix(size(Z_ext,1)/2)+1+fy_range,:)=1;

Z_ext = SurfaceExtension_GerchbergPapoulis(Z_ext, G, Gy, Gox, Goy);

end