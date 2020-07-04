clc;
close all;
clear;

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

load([data_dir 'example_surf_rf.mat']);
pixel_m = median(diff(X(1,:)));
[X, Y, Z] = CleanNaNFromSurfData(X,Y,Z);

%% Extension area
% Sampling intervals
surf_mpp = median(diff(X(1,:)));    ... surface sampling interval [m/pxl]
tif_mpp = median(diff(X_brf(1,:)));

m = size(Z,1);  ... CA height [pixel]
n = size(Z,2);  ... CA width [pixel]

m_ext = round(tif_mpp*(size(Z_tif, 1))*0.5/surf_mpp);   ... extension size in y [pixel]
n_ext = round(tif_mpp*(size(Z_tif, 2))*0.5/surf_mpp);   ... extension size in x [pixel] 

ca_range.y_s = m_ext + 1;   ca_range.y_e = ca_range.y_s + m - 1;   ... y start & end ids of CA in FA [pixel]
ca_range.x_s = n_ext + 1;   ca_range.x_e = ca_range.x_s + n - 1;   ... x start & end ids of CA in FA [pixel]

[X_ext, Y_ext] = meshgrid(-n_ext:n-1+n_ext, -m_ext:m-1+m_ext);  ...extension grid
X_ext = X_ext * surf_mpp + X(1,1);  ... adjust X grid add X(1,1)
Y_ext = Y_ext * surf_mpp + Y(1,1);  ... adjust Y grid add Y(1,1)

Z_ext = NaN(size(X_ext));   ... mark the Z_ext to NaN
Z_ext(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e) = Z;... fill in the valid data point

%% Fit the edge values
% Z_ext(:, 1) = 0;
% Z_ext(1, :) = 0;
% Z_ext(:, end) = 0;
% Z_ext(end, :) = 0;

%% Obtain polynomial bases
order = 64;
X_nor = -1 + 2.*(X_ext - min(X_ext(:)))./(max(X_ext(:)) - min(X_ext(:)));
Y_nor = -1 + 2.*(Y_ext - min(Y_ext(:)))./(max(Y_ext(:)) - min(Y_ext(:)));


[~, ~, ~, z3, ~, ~, ~, ~] = ChebyshevXYJC(X_nor, Y_nor, 1:order, ones(order,1));
% [~, ~, ~, z3, ~, ~, ~, ~] = LegendreXYJC(X_ext, Y_ext, 1:order, ones(order,1));
z3 = reshape(z3, [],size(z3,3));

A = z3(~isnan(Z_ext(:)),:);
b = Z_ext(~isnan(Z_ext(:)));

c = A\b;

[Z_fit,~, ~, ~, ~, ~, ~, ~] = ChebyshevXYJC(X_nor, Y_nor, 1:order, c);

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


