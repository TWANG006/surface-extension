clc;
close all;
clear;

addpath(genpath('../lib/'));

data_dir = '../data/';
brf_dir = '../data/';

%% Load Initial surf & brf
% load brf data
load([brf_dir 'step_0_fluid_jet_tif.mat']);
X_brf = X;
Y_brf = Y;
Z_avg = Z;
brf_params.A = 125e-9/20;
brf_params.sigma_xy = FWHM2Sigma([4.293e-3, 4.293e-3]);
brf_params.d_pix = size(Z_avg, 1);
brf_params.d = brf_params.d_pix * m_per_pixel;
brf_params.lat_res_brf = m_per_pixel;

load([data_dir 'example_surf_rf.mat']);
pixel_m = median(diff(X(1,:)));

% Clean the data
[X, Y, Z] = CleanNaNFromSurfData(X,Y,Z);

%% Dwell time calculation parameters
dx_ibf = 1e-3;
brf_mode = 'avg';
ratio = 1;
tmin=0;
tmax=1;
run_iter = true;
TT = 0;

options = struct(...
    'Algorithm', 'Iterative-FFT',...
    'maxIters', 20, ...
    'PV_dif', 0.01e-9, ...[m]
    'RMS_dif', 0.02e-9, ...[m]
    'dwellTime_dif', 60, ...[s]
    'isDownSampling', false, ...
    'samplingInterval', dx_ibf ... [m]
);


%% 1. Zero extension
[X_ext, Y_ext, Z_0, ca_range] = Surface_Extension(X,Y,Z,brf_params,Z_avg,'zero',false);

Z_ext = NaN(size(X_ext));   ... mark the Z_ext to NaN
Z_ext(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e) = Z;... fill in the valid data point
r = max(size(Z_ext)-size(Z))/2;
% obtain the are of extension
[u,v] = meshgrid(-r:r,-r:r);
[~, rr] = cart2pol(u,v);
se = rr<=r;
BW_Z = imdilate(~isnan(Z_ext), se);
id_ext = BW_Z ==0;

[B... BRF
    , X_B, Y_B...BRF Coordinates
    , ~, ~... full aperture results [m]
    , T_0... dwell time on the dwell grid [s]
    , ~...
    , X_P, Y_P...dwell grid 
    , ~, ~, ~ ... dwell grid coordinates [m]
    , ~, ~, ~...dwell grid results [m]
    , X_ca, Y_ca... clear aperture coordinates [m]
    , ~, ~, ~...[m]
    ] = DwellTime2D_FFT_Full_Test...
    ( X_ext, Y_ext, Z_0... height to remove [m]
    , 0 ...
    , brf_params... BRF parameters
    , brf_mode...
    , X_brf, Y_brf, Z_avg...
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , tmin,tmax...
    , options...
    , ratio...
    , false...
    );

T_0(id_ext) = 0;
Z_removal_dw = ConvFFT2(T_0,B);
Z_0 = RemoveSurface1(X_ext, Y_ext, Z_0);
Z_0 = Z_0 - nanmin(Z_0(:));
Z_residual_dw = Z_0 - nanmin(Z_0(:)) - Z_removal_dw;
Z_residual_ca_0 = Z_residual_dw(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_residual_ca_0 = RemoveSurface1(X_ca, Y_ca, Z_residual_ca_0);

%% 2. Gaussian extension
[~, ~, Z_gauss, ~] = Surface_Extension(X,Y,Z,brf_params,Z_avg,'gauss',false);

[~... BRF
    , ~, ~...BRF Coordinates
    , ~, ~... full aperture results [m]
    , T_gauss... dwell time on the dwell grid [s]
    , ~...
    , ~, ~...dwell grid 
    , ~, ~, ~ ... dwell grid coordinates [m]
    , ~, ~, ~...dwell grid results [m]
    , ~, ~... clear aperture coordinates [m]
    , ~, ~, ~...[m]
    ] = DwellTime2D_FFT_Full_Test...
    ( X_ext, Y_ext, Z_gauss... height to remove [m]
    , 0 ...
    , brf_params... BRF parameters
    , brf_mode...
    , X_brf, Y_brf, Z_avg...
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , tmin,tmax...
    , options...
    , ratio...
    , false...
    );

T_gauss(id_ext) = 0;
Z_removal_dw = ConvFFT2(T_gauss,B);
Z_gauss = RemoveSurface1(X_ext, Y_ext, Z_gauss);
Z_gauss = Z_gauss - nanmin(Z_gauss(:));
Z_residual_dw = Z_gauss - Z_removal_dw;
Z_residual_ca_gauss = Z_residual_dw(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_residual_ca_gauss = RemoveSurface1(X_ca, Y_ca, Z_residual_ca_gauss);

%% 3. 8NN extension
[~, ~, Z_8nn, ~] = Surface_Extension(X,Y,Z,brf_params,Z_avg,'8nn',false);

[~... BRF
    , ~, ~...BRF Coordinates
    , ~, ~... full aperture results [m]
    , T_8nn... dwell time on the dwell grid [s]
    , ~...
    , ~, ~...dwell grid 
    , ~, ~, ~ ... dwell grid coordinates [m]
    , ~, ~, ~...dwell grid results [m]
    , ~, ~... clear aperture coordinates [m]
    , ~, ~, ~...[m]
    ] = DwellTime2D_FFT_Full_Test...
    ( X_ext, Y_ext, Z_8nn... height to remove [m]
    , 0 ...
    , brf_params... BRF parameters
    , brf_mode...
    , X_brf, Y_brf, Z_avg...
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , tmin,tmax...
    , options...
    , ratio...
    , false...
    );

T_8nn(id_ext) = 0;
Z_removal_dw = ConvFFT2(T_8nn,B);
Z_8nn = RemoveSurface1(X_ext, Y_ext, Z_8nn);
Z_8nn = Z_8nn - nanmin(Z_8nn(:));
Z_residual_dw = Z_8nn - Z_removal_dw;
Z_residual_ca_8nn = Z_residual_dw(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_residual_ca_8nn = RemoveSurface1(X_ca, Y_ca, Z_residual_ca_8nn);

%% 4. 8NN extension with fall
[~, ~, Z_8nn_fall, ~] = Surface_Extension(X,Y,Z,brf_params,Z_avg,'8nn',true);

[~... BRF
    , ~, ~...BRF Coordinates
    , ~, ~... full aperture results [m]
    , T_8nn_fall... dwell time on the dwell grid [s]
    , ~...
    , ~, ~...dwell grid 
    , ~, ~, ~ ... dwell grid coordinates [m]
    , ~, ~, ~...dwell grid results [m]
    , ~, ~... clear aperture coordinates [m]
    , ~, ~, ~...[m]
    ] = DwellTime2D_FFT_Full_Test...
    ( X_ext, Y_ext, Z_8nn_fall... height to remove [m]
    , 0 ...
    , brf_params... BRF parameters
    , brf_mode...
    , X_brf, Y_brf, Z_avg...
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , tmin,tmax...
    , options...
    , ratio...
    , false...
    );

T_8nn_fall(id_ext) = 0;
Z_removal_dw = ConvFFT2(T_8nn_fall,B);
Z_8nn_fall = RemoveSurface1(X_ext, Y_ext, Z_8nn_fall);
Z_8nn_fall = Z_8nn_fall - nanmin(Z_8nn_fall(:));
Z_residual_dw = Z_8nn_fall - Z_removal_dw;
Z_residual_ca_8nn_fall = Z_residual_dw(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_residual_ca_8nn_fall = RemoveSurface1(X_ca, Y_ca, Z_residual_ca_8nn_fall);

%% 5. Smooth extension
[~, ~, Z_smooth, ~] = Surface_Extension(X,Y,Z,brf_params,Z_avg,'smooth',false);

[~... BRF
    , ~, ~...BRF Coordinates
    , ~, ~... full aperture results [m]
    , T_smooth... dwell time on the dwell grid [s]
    , ~...
    , ~, ~...dwell grid 
    , ~, ~, ~ ... dwell grid coordinates [m]
    , ~, Z_removal_dw_smooth, ~...dwell grid results [m]
    , ~, ~... clear aperture coordinates [m]
    , ~, ~, ~...[m]
    ] = DwellTime2D_FFT_Full_Test...
    ( X_ext, Y_ext, Z_smooth... height to remove [m]
    , 0 ...
    , brf_params... BRF parameters
    , brf_mode...
    , X_brf, Y_brf, Z_avg...
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , tmin,tmax...
    , options...
    , ratio...
    , false...
    );

T_smooth(id_ext) = 0;
Z_removal_dw = ConvFFT2(T_smooth,B);
Z_smooth = RemoveSurface1(X_ext, Y_ext, Z_smooth);
Z_smooth = Z_smooth - nanmin(Z_smooth(:));
Z_residual_dw = Z_smooth - Z_removal_dw;
Z_residual_ca_smooth = Z_residual_dw(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_residual_ca_smooth = RemoveSurface1(X_ca, Y_ca, Z_residual_ca_smooth);


%% 6. Smooth extension with fall
[~, ~, Z_smooth_fall, ~] = Surface_Extension(X,Y,Z,brf_params,Z_avg,'smooth',true);

[~... BRF
    , ~, ~...BRF Coordinates
    , ~, ~... full aperture results [m]
    , T_smooth_fall... dwell time on the dwell grid [s]
    , ~...
    , ~, ~...dwell grid 
    , ~, ~, ~ ... dwell grid coordinates [m]
    , ~, ~, ~...dwell grid results [m]
    , ~, ~... clear aperture coordinates [m]
    , ~, ~, ~...[m]
    ] = DwellTime2D_FFT_Full_Test...
    ( X_ext, Y_ext, Z_smooth_fall... height to remove [m]
    , 0 ...
    , brf_params... BRF parameters
    , brf_mode...
    , X_brf, Y_brf, Z_avg...
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , tmin,tmax...
    , options...
    , ratio...
    , false...
    );

T_smooth_fall(id_ext) = 0;
Z_removal_dw = ConvFFT2(T_smooth_fall,B);
Z_smooth_fall = RemoveSurface1(X_ext, Y_ext, Z_smooth_fall);
Z_smooth_fall = Z_smooth_fall - nanmin(Z_smooth_fall(:));
Z_residual_dw = Z_smooth_fall - Z_removal_dw;
Z_residual_ca_smooth_fall = Z_residual_dw(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_residual_ca_smooth_fall = RemoveSurface1(X_ca, Y_ca, Z_residual_ca_smooth_fall);

%% 7. Gerchberg extension
[~, ~, Z_gp, ~] = Surface_Extension(X,Y,Z,brf_params,Z_avg,'gerchberg',false,-18:18,-3:3);

[~... BRF
    , ~, ~...BRF Coordinates
    , ~, ~... full aperture results [m]
    , T_gp... dwell time on the dwell grid [s]
    , ~...
    , ~, ~...dwell grid 
    , ~, ~, ~ ... dwell grid coordinates [m]
    , ~, ~, ~...dwell grid results [m]
    , ~, ~... clear aperture coordinates [m]
    , ~, ~, ~...[m]
    ] = DwellTime2D_FFT_Full_Test...
    ( X_ext, Y_ext, Z_gp... height to remove [m]
    , 0 ...
    , brf_params... BRF parameters
    , brf_mode...
    , X_brf, Y_brf, Z_avg...
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , tmin,tmax...
    , options...
    , ratio...
    , false...
    );

T_gp(id_ext) = 0;
Z_removal_dw = ConvFFT2(T_gp,B);
Z_gp = RemoveSurface1(X_ext, Y_ext, Z_gp);
Z_gp = Z_gp - nanmin(Z_gp(:));
Z_residual_dw = Z_gp - Z_removal_dw;
Z_removal_ca = Z_removal_dw(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_residual_ca_gp = Z_residual_dw(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_residual_ca_gp = RemoveSurface1(X_ca, Y_ca, Z_residual_ca_gp);

%% 8. Run iterative extension
if(run_iter)
    [T_iter, ~, Z_residual_ca_iter] = Surface_Extension_Iter(X,Y,Z,brf_params,'avg',1e-3, 0.7e-9,X_brf,Y_brf,Z_avg,...
        'smooth',false,[],[],[],[],[],...
        'smooth',false,[],[],[],[],[]);
end

%% Show the results
fsfig('One-step surface extension dwell time results');
subplot(9,3,3);
surf(X*1e3,Y*1e3,Z*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
rms_Z = std(Z(:),1)*1e9;
title(['Original Surface: PV =' num2str(round((max(Z(:))-min(Z(:)))*1e9,2)) ' nm, RMS = ' num2str(round(rms_Z,2)) ' nm']);
axis image

subplot(9,3,4)
surf(X_ext*1e3,Y_ext*1e3,Z_0*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
title('Zero extension');
axis image
subplot(9,3,5)
surf(X_ext*1e3,Y_ext*1e3,T_0,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[s]';
title(['Zero extension: Total dwell time = ' num2str(round(sum(T_0(:))/3600),2) ' h']);
axis image
subplot(9,3,6)
surf(X*1e3,Y*1e3,Z_residual_ca_0*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
rms_Z = std(Z_residual_ca_0(:),1)*1e9;
caxis([-1 1]*3*rms_Z);
title(['Residual: PV =' num2str(round((max(Z_residual_ca_0(:))-min(Z_residual_ca_0(:)))*1e9,2)) ' nm, RMS = ' num2str(round(rms_Z,2)) ' nm']);
axis image


subplot(9,3,7)
surf(X_ext*1e3,Y_ext*1e3,Z_gauss*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
title('Gaussian extension');
axis image
subplot(9,3,8)
surf(X_ext*1e3,Y_ext*1e3,T_gauss,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[s]';
title(['Gaussian extension: Total dwell time = ' num2str(round(sum(T_gauss(:))/3600),2) ' h']);
axis image
subplot(9,3,9)
surf(X*1e3,Y*1e3,Z_residual_ca_gauss*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
rms_Z = std(Z_residual_ca_gauss(:),1)*1e9;
caxis([-1 1]*3*rms_Z);
title(['Residual: PV =' num2str(round((max(Z_residual_ca_gauss(:))-min(Z_residual_ca_gauss(:)))*1e9,2)) ' nm, RMS = ' num2str(round(rms_Z,2)) ' nm']);
axis image

subplot(9,3,10)
surf(X_ext*1e3,Y_ext*1e3,Z_8nn*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
title('8NN extension');
axis image
subplot(9,3,11)
surf(X_ext*1e3,Y_ext*1e3,T_8nn,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[s]';
title(['8NN extension: Total dwell time = ' num2str(round(sum(T_8nn(:))/3600),2) ' h']);
axis image
subplot(9,3,12)
surf(X*1e3,Y*1e3,Z_residual_ca_8nn*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
rms_Z = std(Z_residual_ca_8nn(:),1)*1e9;
caxis([-1 1]*3*rms_Z);
title(['Residual: PV =' num2str(round((max(Z_residual_ca_8nn(:))-min(Z_residual_ca_8nn(:)))*1e9,2)) ' nm, RMS = ' num2str(round(rms_Z,2)) ' nm']);
axis image

subplot(9,3,13)
surf(X_ext*1e3,Y_ext*1e3,Z_8nn_fall*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
title('8NN extension with fall');
axis image
subplot(9,3,14)
surf(X_ext*1e3,Y_ext*1e3,T_8nn_fall,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[s]';
title(['8NN extension with fall: Total dwell time = ' num2str(round(sum(T_8nn_fall(:))/3600),2) ' h']);
axis image
subplot(9,3,15)
surf(X*1e3,Y*1e3,Z_residual_ca_8nn_fall*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
rms_Z = std(Z_residual_ca_8nn_fall(:),1)*1e9;
caxis([-1 1]*3*rms_Z);
title(['Residual: PV =' num2str(round((max(Z_residual_ca_8nn_fall(:))-min(Z_residual_ca_8nn_fall(:)))*1e9,2)) ' nm, RMS = ' num2str(round(rms_Z,2)) ' nm']);
axis image

subplot(9,3,16)
surf(X_ext*1e3,Y_ext*1e3,Z_smooth*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
title('Smooth extension');
axis image
subplot(9,3,17)
surf(X_ext*1e3,Y_ext*1e3,T_smooth,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[s]';
title(['Smooth extension: Total dwell time = ' num2str(round(sum(T_smooth(:))/3600),2) ' h']);
axis image
subplot(9,3,18)
surf(X*1e3,Y*1e3,Z_residual_ca_smooth*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
rms_Z = std(Z_residual_ca_smooth(:),1)*1e9;
caxis([-1 1]*3*rms_Z);
title(['Residual: PV =' num2str(round((max(Z_residual_ca_smooth(:))-min(Z_residual_ca_smooth(:)))*1e9,2)) ' nm, RMS = ' num2str(round(rms_Z,2)) ' nm']);
axis image

subplot(9,3,19)
surf(X_ext*1e3,Y_ext*1e3,Z_smooth_fall*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
title('Smooth extension with fall');
axis image
subplot(9,3,20)
surf(X_ext*1e3,Y_ext*1e3,T_smooth_fall,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[s]';
title(['Smooth extension with fall: Total dwell time = ' num2str(round(sum(T_smooth_fall(:))/3600),2) ' h']);
axis image
subplot(9,3,21)
surf(X*1e3,Y*1e3,Z_residual_ca_smooth_fall*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
rms_Z = std(Z_residual_ca_smooth_fall(:),1)*1e9;
caxis([-1 1]*3*rms_Z);
title(['Residual: PV =' num2str(round((max(Z_residual_ca_smooth_fall(:))-min(Z_residual_ca_smooth_fall(:)))*1e9,2)) ' nm, RMS = ' num2str(round(rms_Z,2)) ' nm']);
axis image

subplot(9,3,22)
surf(X_ext*1e3,Y_ext*1e3,Z_gp*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
title('Gerchberg extension');
axis image
subplot(9,3,23)
surf(X_ext*1e3,Y_ext*1e3,T_gp,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[s]';
title(['Gerchberg extension: Total dwell time = ' num2str(round(sum(T_gp(:))/3600),2) ' h']);
axis image
subplot(9,3,24)
surf(X*1e3,Y*1e3,Z_residual_ca_gp*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
rms_Z = std(Z_residual_ca_gp(:),1)*1e9;
caxis([-1 1]*3*rms_Z);
title(['Residual: PV =' num2str(round((max(Z_residual_ca_gp(:))-min(Z_residual_ca_gp(:)))*1e9,2)) ' nm, RMS = ' num2str(round(rms_Z,2)) ' nm']);
axis image

subplot(9,3,26)
surf(X_ext*1e3,Y_ext*1e3,T_iter,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[s]';
title(['Iterative extension: Total dwell time = ' num2str(round(sum(T_iter(:))/3600),2) ' h']);
axis image
subplot(9,3,27)
surf(X*1e3,Y*1e3,Z_residual_ca_iter*1e9,'EdgeColor', 'none');view([0 90]); 
c = colorbar;
c.Label.String = '[nm]';
rms_Z = std(Z_residual_ca_iter(:),1)*1e9;
caxis([-1 1]*3*rms_Z);
title(['Residual: PV =' num2str(round((max(Z_residual_ca_iter(:))-min(Z_residual_ca_iter(:)))*1e9,2)) ' nm, RMS = ' num2str(round(rms_Z,2)) ' nm']);
axis image