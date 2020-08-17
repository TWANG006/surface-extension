function [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension_Smooth(...
    X, Y, Z,...unextended surface error map
    tif_mpp,...TIF sampling interval [m/pxl]
    Z_tif...TIF profile
    )
% Function
%   [X_ext, Y_ext, Z_ext] = Surface_Extenstion_Spline(X, Y, Z, X_tif, Y_tif, Z_tif)
% Purpose
%   Extend the surface error map using spline interpolation

%% 0. Obtain required parameters
% Sampling intervals
surf_mpp = median(diff(X(1,:)));    ... surface sampling interval [m/pxl]

m = size(Z,1);  ... CA height [pixel]
n = size(Z,2);  ... CA width [pixel]

m_ext = round(tif_mpp*(size(Z_tif, 1))*0.5/surf_mpp);   ... extension size in y [pixel]
n_ext = round(tif_mpp*(size(Z_tif, 2))*0.5/surf_mpp);   ... extension size in x [pixel] 
r = max(m_ext, n_ext);

ca_range.v_s = m_ext + 1;   ca_range.v_e = ca_range.v_s + m - 1;   ... y start & end ids of CA in FA [pixel]
ca_range.u_s = n_ext + 1;   ca_range.u_e = ca_range.u_s + n - 1;   ... x start & end ids of CA in FA [pixel]


%% 1. Initial extension matrices
% extension sizes
[X_ext, Y_ext] = meshgrid(-n_ext:n-1+n_ext, -m_ext:m-1+m_ext);  ...extension grid
X_ext = X_ext * surf_mpp + X(1,1);  ... adjust X grid add X(1,1)
Y_ext = Y_ext * surf_mpp + Y(1,1);  ... adjust Y grid add Y(1,1)
Z_ext = NaN(size(X_ext));   ... mark the Z_ext to NaN
Z_ext(ca_range.v_s:ca_range.v_e, ca_range.u_s:ca_range.u_e) = Z;... fill in the valid data point
id = ~isnan(Z);
% id_invalid  = isnan(Z_ext);

% obtain the are of extension
[u,v] = meshgrid(-r:r,-r:r);
[~, rr] = cart2pol(u,v);
se = rr<=r;
BW_Z = imdilate(~isnan(Z_ext), se);
id_ext = BW_Z ==1;


% u_edge=[X_ext(1, 1);X_ext(1,end)';X_ext(end,1);X_ext(end,end)']; 
% v_edge=[Y_ext(1, 1);Y_ext(1,end)';Y_ext(end, 1);Y_ext(end,end)'];
% z_edge = 0*u_edge;
% p = 8;
% u_edge=[X_ext(1:p:end, 1);X_ext(1,1:p:end)';X_ext(1:p:end,end);X_ext(end,1:p:end)']; 
% v_edge=[Y_ext(1:p:end, 1);Y_ext(1,1:p:end)';Y_ext(1:p:end,end);Y_ext(end,1:p:end)'];
% z_edge = 0*u_edge;


F = scatteredInterpolant(X(id),Y(id),Z(id), 'natural');
% Z_ext(id_ext) = F(X_ext(id_ext), Y_ext(id_ext));


% Z_ext(isnan(Z_ext)) = 0;
Z_ext = F(X_ext, Y_ext);

% Z_ext(id_invalid) = F(X_ext(id_invalid), Y_ext(id_invalid));
% Z_ext = F(X_ext, Y_ext);

% Z_ext = NaN(size(X_ext));   ... mark the Z_ext to NaN
% Z_ext(ca_range.v_s:ca_range.v_e, ca_range.u_s:ca_range.u_e) = Z;... fill in the valid data points
% BW_ini = ~isnan(Z_ext);   ... obtain the black&white map
% BW_prev = BW_ini;
% % BW_all=0;
% h = size(Z_ext, 1);
% w = size(Z_ext, 2);



end