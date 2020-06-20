clc;
close all;
clear;

addpath(genpath('../lib/'));

data_dir = '../data/';
brf_dir = '../data/';

%% Load Initial surf & brf (beam removal function)
% load brf data
load([brf_dir 'step_0_fluid_jet_tif.mat']);
X_tif = X;
Y_tif = Y;
Z_tif = Z;
brf_params.A = 125e-9/20;
brf_params.sigma_xy = FWHM2Sigma([4.293e-3, 4.293e-3]);
brf_params.d_pix = size(Z_tif, 1);
brf_params.d = brf_params.d_pix * m_per_pixel;
brf_params.lat_res_brf = m_per_pixel;
pixel_m = median(diff(X(1,:)));

% load([data_dir 'example_surf_cf.mat']);
load([data_dir 'example_surf_rf.mat']);


%% Call different surface extesion algorithm
[X_ext, Y_ext, Z_smth, ca_range] = Surface_Extension(X,Y,Z,brf_params,Z_tif,'smooth', false);
[~, ~, Z_smth_fall, ~] = Surface_Extension(X,Y,Z,brf_params,Z_tif, 'smooth', true);
[~, ~, Z_8nn, ~] = Surface_Extension(X,Y,Z,brf_params,Z_tif, '8nn', false);
[~, ~, Z_8nn_fall, ~] = Surface_Extension(X,Y,Z,brf_params,Z_tif, '8nn', true);
[~, ~, Z_gauss, ~] = Surface_Extension(X,Y,Z,brf_params,Z_tif,'gauss',false);
[~, ~, Z_gp, ~] = Surface_Extension(X,Y,Z,brf_params,Z_tif,'gerchberg',false,-18:18,-3:3);
[~, ~, Z_0, ~] = Surface_Extension(X,Y,Z,brf_params,Z_tif,'zero',false);

x = X(round(size(X,1)*0.5),:);
z = Z(round(size(X,1)*0.5),:);

x_ext = X_ext(round(size(X_ext,1)*0.5),:);
z_smth = Z_smth(round(size(X_ext,1)*0.5),:);
z_smth_fall = Z_smth_fall(round(size(X_ext,1)*0.5),:);
z_8nn = Z_8nn(round(size(X_ext,1)*0.5),:);
z_8nn_fall = Z_8nn_fall(round(size(X_ext,1)*0.5),:);
z_gauss = Z_gauss(round(size(X_ext,1)*0.5),:);
z_gp = Z_gp(round(size(X_ext,1)*0.5),:);
z_0 = Z_0(round(size(X_ext,1)*0.5),:);


% cut the edges of the area of extension
Z_ext = NaN(size(X_ext));   ... mark the Z_ext to NaN
Z_ext(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e) = Z;... fill in the valid data point
r = max(size(Z_ext)-size(Z))/2;

[u,v] = meshgrid(-r:r,-r:r);
[~, rr] = cart2pol(u,v);
se = rr<=r;
BW_Z = imdilate(~isnan(Z_ext), se);
id_ext = BW_Z ==0;

Z_8nn(id_ext) = NaN;
Z_8nn_fall(id_ext) = NaN;
Z_smth(id_ext) = NaN;
Z_smth_fall(id_ext) = NaN;
Z_gauss(id_ext) = NaN;
Z_gp(id_ext) = NaN;
Z_0(id_ext) = NaN;

%% Show the results
fsfig('');
subplot(7,2,1)
surf(X_ext*1e3,Y_ext*1e3,Z_0*1e9,'EdgeColor', 'none');view([0 90]);axis image; colorbar;
title('Zero extension');
subplot(7,2,3)
surf(X_ext*1e3,Y_ext*1e3,Z_gauss*1e9,'EdgeColor', 'none');view([0 90]);axis image; colorbar;
title('Gaussian extension');
subplot(7,2,5)
surf(X_ext*1e3,Y_ext*1e3,Z_8nn*1e9,'EdgeColor', 'none');view([0 90]);axis image;colorbar;
title('8-Nearest-Neighbors extension');
subplot(7,2,7)
surf(X_ext*1e3,Y_ext*1e3,Z_8nn_fall*1e9,'EdgeColor', 'none');view([0 90]);axis image;colorbar;
title('8-Nearest-Neighbors extension + fall');
subplot(7,2,9)
surf(X_ext*1e3,Y_ext*1e3,Z_smth*1e9,'EdgeColor', 'none');view([0 90]);axis image;colorbar;
title('Smooth extension');
subplot(7,2,11)
surf(X_ext*1e3,Y_ext*1e3,Z_smth_fall*1e9,'EdgeColor', 'none');view([0 90]);axis image;colorbar;
title('Smooth extension + fall');
subplot(7,2,13)
surf(X_ext*1e3,Y_ext*1e3,Z_gp*1e9,'EdgeColor', 'none');view([0 90]);axis image;colorbar;
title('Gerchberg Pouplis extension');

subplot(7,2,4)
surf(X*1e3,Y*1e3,Z*1e9,'EdgeColor', 'none');view([0 90]);axis image;colorbar;
title('Original Surface');

Colors = lines(7);
subplot(7,2,[8,10,12]);
plot(x_ext, z_smth, 'LineWidth', 1);hold on;
plot(x_ext, z_smth_fall, 'LineWidth', 1);hold on;
plot(x_ext, z_8nn, 'LineWidth', 1);hold on;
plot(x_ext, z_8nn_fall, 'LineWidth', 1);hold on;
plot(x_ext, z_gauss, 'LineWidth', 1); hold on;
plot(x_ext, z_gp, 'LineWidth', 1); hold on;
plot(x_ext, z_0, 'LineWidth', 1); hold on;
plot(x, z, 'LineWidth', 2);hold off;

% plot(x_ext, z_smth, 'Color', Colors(1,:),'LineWidth', 2);hold on;
% plot(x_ext, z_smth_fall, 'Color', Colors(2,:),'LineWidth', 2);hold on;
% plot(x_ext, z_8nn, 'Color', Colors(3,:),'LineWidth', 2);hold on;
% plot(x_ext, z_8nn_fall, 'Color', Colors(4,:),'LineWidth', 2);hold on;
% plot(x_ext, z_gauss, 'Color', Colors(5,:), 'LineWidth', 2); hold on;
% plot(x_ext, z_gp, 'Color', Colors(6,:), 'LineWidth', 2); hold on;
% plot(x, z, 'Color', Colors(7,:),'LineWidth', 2);hold off;

legend('smooth', 'smooth fall', 'flat', 'flat fall', 'Gauss', 'Gerchberg', 'original');
title('Profiles of the center line');
