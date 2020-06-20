function [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension_Gauss(...
    X, Y, Z,...unextended surface error map
    brf_params,...brf parameters [m/pxl]
    Z_tif...TIF profile
    )
% Function
%   [X_ext, Y_ext, Z_ext] = Surface_Extenstion_Spline(X, Y, Z, X_tif, Y_tif, Z_tif)
% Purpose
%   Extend the surface error map using spline interpolation

%% 0. Obtain required parameters
% Sampling intervals
surf_mpp = median(diff(X(1,:)));    ... surface sampling interval [m/pxl]
sigma = brf_params.sigma_xy(1);

m = size(Z,1);  ... CA height [pixel]
n = size(Z,2);  ... CA width [pixel]

m_ext = round(brf_params.lat_res_brf*(size(Z_tif, 1))*0.5/surf_mpp);   ... extension size in y [pixel]
n_ext = round(brf_params.lat_res_brf*(size(Z_tif, 2))*0.5/surf_mpp);   ... extension size in x [pixel] 

ca_range.y_s = m_ext + 1;   ca_range.y_e = ca_range.y_s + m - 1;   ... y start & end ids of CA in FA [pixel]
ca_range.x_s = n_ext + 1;   ca_range.x_e = ca_range.x_s + n - 1;   ... x start & end ids of CA in FA [pixel]


%% 1. Initial extension matrices
% extension sizes
[X_ext, Y_ext] = meshgrid(-n_ext:n-1+n_ext, -m_ext:m-1+m_ext);  ...extension grid
X_ext = X_ext * surf_mpp + X(1,1);  ... adjust X grid add X(1,1)
Y_ext = Y_ext * surf_mpp + Y(1,1);  ... adjust Y grid add Y(1,1)
Z_ext = NaN(size(X_ext));   ... mark the Z_ext to NaN
Z_ext(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e) = Z;... fill in the valid data point


%% Finding edge points
id_edg = Surface_Extension_EdgeExtraction(Z_ext);
x_edg = X_ext(id_edg);
y_edg = Y_ext(id_edg);
z_edg = Z_ext(id_edg);

id_fil = isnan(Z_ext); ... filled data ids
x_fil = X_ext(id_fil);  ... x coordinates of filled data
y_fil = Y_ext(id_fil);  ... y coordinates of filled data

%% Calculate the gaussian profile
gauss_profiles = zeros*x_fil;
for k = 1:length(x_fil)
    % min distances from filled points to edge points
    [min_dist, i] = min(sqrt((x_fil(k) - x_edg).^2+(y_fil(k) - y_edg).^2));
    
    % calculate the fall profile
    gauss_profiles(k) = z_edg(i) * exp(-min_dist.^2/(2*sigma.^2));
end

Z_ext(id_fil) = gauss_profiles;


end