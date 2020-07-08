clc;
% close all;
clear;

addpath(genpath('../lib/'));

data_dir = '../data/';
brf_dir = '../data/';

%% Load Initial surf & brf
% load brf data
load([brf_dir 'step_0_fluid_jet_tif.mat']);
X_brf = X;
Y_brf = Y;
Z_tif = Z;
brf_params.A = 125e-9/20;
brf_params.sigma_xy = FWHM2Sigma([4.293e-3, 4.293e-3]);
brf_params.d_pix = size(Z_tif, 1);
brf_params.d = brf_params.d_pix * m_per_pixel;
brf_params.lat_res_brf = m_per_pixel;
% load([brf_dir 'example_brf.mat']);

load([data_dir 'example_surf_rf.mat']);
pixel_m = median(diff(X(1,:)));
[X, Y, Z] = CleanNaNFromSurfData(X,Y,Z);

%% Extension area
% Sampling intervals
% surf_mpp = median(diff(X(1,:)));    ... surface sampling interval [m/pxl]
tif_mpp = median(diff(X_brf(1,:)));
% r_pix = round(tif_mpp*(size(Z_tif, 1))*0.5/surf_mpp);
% 
% m = size(Z,1);  ... CA height [pixel]
% n = size(Z,2);  ... CA width [pixel]
% 
% m_ext = round(tif_mpp*(size(Z_tif, 1))*0.5/surf_mpp);   ... extension size in y [pixel]
% n_ext = round(tif_mpp*(size(Z_tif, 2))*0.5/surf_mpp);   ... extension size in x [pixel] 
% 
% % clear aperture
% ca_range.y_s = m_ext + 1;   ca_range.y_e = ca_range.y_s + m - 1;   ... y start & end ids of CA in FA [pixel]
% ca_range.x_s = n_ext + 1;   ca_range.x_e = ca_range.x_s + n - 1;   ... x start & end ids of CA in FA [pixel]
% 
% % dwell grid
% dw_range.y_s = ca_range.y_s - r_pix;
% dw_range.x_s = ca_range.x_s - r_pix;
% dw_range.y_e = ca_range.y_e + r_pix;
% dw_range.x_e = ca_range.x_e + r_pix;
% 
% 
% [X_ext, Y_ext] = meshgrid(-n_ext:n-1+n_ext, -m_ext:m-1+m_ext);  ...extension grid
% X_ext = X_ext * surf_mpp + X(1,1);  ... adjust X grid add X(1,1)
% Y_ext = Y_ext * surf_mpp + Y(1,1);  ... adjust Y grid add Y(1,1)
% 
% Z_ext = NaN(size(X_ext));   ... mark the Z_ext to NaN
% Z_ext(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e) = Z;... fill in the valid data point

%% Fit the edge values
% zero boundary contidion
% w = 10;
% Z_ext(:, 1) = 0;
% Z_ext(1, :) = 0;
% Z_ext(:, end) = 0;
% Z_ext(end, :) = 0;
% 
% W = ones(size(Z_ext));
% W(:, 1) = w;
% W(1, :) = w;
% W(:, end) = w;
% W(end, :) = w;
% % smooth boundary condition
% [X_ext, Y_ext, Z_smth, ca_range] = Surface_Extension(X,Y,Z,brf_params,Z_tif,'smooth', false);
% Z_ext(:, 1) = Z_smth(:, 1);
% Z_ext(1, :) = Z_smth(1, :);
% Z_ext(:, end) = Z_smth(:, end);
% Z_ext(end, :) = Z_smth(end, :);

%% Obtain polynomial fitting
[X_ext, Y_ext, Z_fit, ca_range] = Surface_Extension(...
    X, Y, Z,...unextended surface error map
    brf_params,...TIF sampling interval [m/pxl]
    Z_tif,...TIF profile
    'poly',...
    false,...
    [],[],...
    8, 256,...polynomial orders in y, x
    'Chebyshev'...Chebyshev or Legendre
    );

%% Show fitting result
fsfig('');
subplot(3,2,[1,3,5]);
surf(X_ext*1e3, Y_ext*1e3, 1e9*Z_fit, 'EdgeColor', 'none');
view([0 90]);
axis image;
colorbar;

subplot(322);
surf(X*1e3, Y*1e3, 1e9*Z_fit(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e), 'EdgeColor', 'none');
view([0 90]);
axis image;
colorbar;

subplot(324);
surf(X*1e3, Y*1e3, 1e9*Z, 'EdgeColor', 'none');
view([0 90]);
axis image;
colorbar;

subplot(326);
Z_residual = 1e9*(Z_fit(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e) - Z);
surf(X*1e3, Y*1e3, Z_residual, 'EdgeColor', 'none');
view([0 90]);
colorbar;
axis image;
title(['Residual = ' num2str(nanstd(Z_residual(:),1)) ' nm']);

%% calculation
dx_ibf = 1e-3;
brf_mode = 'avg';
ratio = 1;
tmin=0;
tmax=1;
run_iter = true;
TT = 0;

options = struct(...
    'Algorithm', 'Iterative-FFT',...
    'maxIters', 50, ...
    'PV_dif', 0.01e-9, ...[m]
    'RMS_dif', 0.02e-9, ...[m]
    'dwellTime_dif', 60, ...[s]
    'isDownSampling', false, ...
    'samplingInterval', dx_ibf ... [m]
);

Z_ext = NaN(size(X_ext));   ... mark the Z_ext to NaN
Z_ext(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e) = Z_fit(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);... fill in the valid data point
r = max(size(Z_ext)-size(Z))/2;
% obtain the are of extension
[u,v] = meshgrid(-r:r,-r:r);
[~, rr] = cart2pol(u,v);
se = rr<=r;
BW_Z = imdilate(~isnan(Z_ext), se);
id_ext = BW_Z ==0;

Z_fit(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e) = Z;
[B... BRF
    , X_B, Y_B...BRF Coordinates
    , ~, ~... full aperture results [m]
    , T... dwell time on the dwell grid [s]
    , ~...
    , X_P, Y_P...dwell grid 
    , ~, ~, ~ ... dwell grid coordinates [m]
    , ~, ~, ~...dwell grid results [m]
    , X_ca, Y_ca... clear aperture coordinates [m]
    , ~, ~, ~...[m]
    ] = DwellTime2D_FFT_Full_Test...
    ( X_ext, Y_ext, Z_fit... height to remove [m]
    , 0 ...
    , brf_params... BRF parameters
    , brf_mode...
    , X_brf, Y_brf, Z_tif...
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , tmin,tmax...
    , options...
    , ratio...
    , false...
    );

T(id_ext) = 0;
Z_removal_dw = ConvFFT2(T,B);
Z_fit_rm = RemoveSurface1(X_ext, Y_ext, Z_fit);
Z_fit_rm = Z_fit_rm - nanmin(Z_fit_rm(:));
Z_residual_dw = Z_fit_rm - nanmin(Z_fit_rm(:)) - Z_removal_dw;
Z_residual_ca = Z_residual_dw(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_residual_ca = RemoveSurface1(X_ca, Y_ca, Z_residual_ca);

%% Iterative solution
[T_iter, B, Z_residual_ca_iter] = Surface_Extension_Iter(...
    X,Y,Z,...
    brf_params,...
    'avg',...
    1e-3, ...
    0.7e-9,...
    X_brf,Y_brf,Z_tif,...
    'poly',false,[],[],8,256,'Chebyshev',...
    'poly',false,[],[],6,250,'Chebyshev');

% [T_iter, B, Z_residual_ca_iter] = Surface_Extension_Iter(...
%     X,Y,Z,...
%     brf_params,...
%     'avg',...
%     1e-3, ...
%     0.7e-9,...
%     X_brf,Y_brf,Z_tif,...
%     'poly',false,[],[],8,256,'Chebyshev',...
%     'smooth',false,[],[],[],[],[]);

%% display result
fsfig('');
subplot(2,2,1)
surf(X_ext*1e3,Y_ext*1e3,T,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[s]';
title(['Polyfitting extension: Total dwell time = ' num2str(round(sum(T(:))/3600),2) ' h']);
axis image
subplot(2,2,2)
surf(X*1e3,Y*1e3,Z_residual_ca*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
rms_Z = std(Z_residual_ca(:),1)*1e9;
caxis([-1 1]*3*rms_Z);
title(['Initial estiated residual: PV =' num2str(round((max(Z_residual_ca(:))-min(Z_residual_ca(:)))*1e9,2)) ' nm, RMS = ' num2str(round(rms_Z,2)) ' nm']);
axis image

subplot(2,2,3)
surf(X_ext*1e3,Y_ext*1e3,T_iter,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[s]';
title(['Iterative extension: Total dwell time = ' num2str(round(sum(T(:))/3600),2) ' h']);
axis image

subplot(2,2,4)
surf(X*1e3,Y*1e3,Z_residual_ca_iter*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
rms_Z = std(Z_residual_ca_iter(:),1)*1e9;
caxis([-1 1]*3*rms_Z);
title(['Residual after iterative extension: PV =' num2str(round((max(Z_residual_ca_iter(:))-min(Z_residual_ca_iter(:)))*1e9,2)) ' nm, RMS = ' num2str(round(rms_Z,2)) ' nm']);
axis image

% Z_removal_dw = ConvFFT2(T_iter, B);
% Z_residual_dw = Z_ext - Z_removal_dw;
% Z_residual_dw = RemoveSurface1(X_ext, Y_ext, Z_residual_dw);
% subplot(3,2,6)
% surf(X_ext*1e3,Y_ext*1e3,Z_residual_dw*1e9,'EdgeColor', 'none');view([0 90]); 
% c = colorbar;
% c.Label.String = '[nm]';
% rms_Z = nanstd(Z_residual_dw(:),1)*1e9;
% caxis([-1 1]*3*rms_Z);
% title(['Residual: PV =' num2str(round((max(Z_residual_dw(:))-min(Z_residual_dw(:)))*1e9,2)) ' nm, RMS = ' num2str(round(rms_Z,2)) ' nm']);
% axis image