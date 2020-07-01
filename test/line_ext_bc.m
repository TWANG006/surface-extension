% line extension using polyfit
% least square fit with boundary condition 
% Heejoo heejoooptics@gmail.com
% 2020.06.29 start..

%% load the surface data
clc;
close all;
clear;

addpath(genpath('../lib/'));

data_dir = '../data/';
brf_dir = '../data/';

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

%% pick up the line and polyfit

line_num = 1; fit_input = Z(:,line_num); fit_inputX = Y(:,line_num);
figure(1);plot(fit_inputX,fit_input);
%% column extension
clear fitResult
for N = 1:size(Z,2)
    line_num = N;
    % scale up
    fit_input = Z(:,line_num).*1e6; fit_inputX = Y(:,line_num).*1e3; dx = fit_inputX(3)-fit_inputX(2);
    % falling/extension boundarycondition
    fit_input = [0 ;fit_input ;0];fit_inputX = [fit_inputX(1)-30*dx; fit_inputX ;fit_inputX(end)+30*dx];
    % n th order polynomial
    Order = 6;
    A = zeros(length(fit_inputX),Order);
    
    for M = 1:Order
    A(:,end-M+1) = fit_inputX.^M;    
    end
    
    A(:,end+1) = ones(size(fit_inputX));
    B = fit_input;
    
    fitParameters = A\B; % least squre fit
    
    
    reConstructionX = fit_inputX(1):fit_inputX(3)-fit_inputX(2):fit_inputX(end);
    fitResult(:,N) = polyval(fitParameters,reConstructionX);
    
%     figure(1);plot(fit_inputX,fit_input,'*');pause(0.1);title('Original Map');
%     figure(2);plot(reConstructionX,fitResult,'*');pause(0.1);title('fitting Map');
%     figure(3);plot(fit_inputX,fit_input-fitResult,'*');pause(0.1);title('residual Map');

end

figure(2);
surf(fitResult*1e3, 'EdgeColor', 'none');
% axis image; 
colormap jet;
c = colorbar;
c.Label.String = '[nm]';
view([0 90]);
axis image;
title('Extended');


figure(4);
subplot(411);
surf(fitResult(31:56,:)*1e3, 'EdgeColor', 'none');
% axis image; 
colormap jet;
c = colorbar;
c.Label.String = '[nm]';
view([0 90]);
axis image;
title('Fitted clear aperture');

subplot(412);
surf(Z*1e9, 'EdgeColor', 'none');
% axis image; 
colormap jet;
c = colorbar;
c.Label.String = '[nm]';
view([0 90]);
axis image;
title('Original clear aperture');

subplot(413);
testR1 = fitResult(31:56,:)*1e3 - Z*1e9;
surf(fitResult(31:56,:)*1e3 - Z*1e9, 'EdgeColor', 'none');
% axis image;
colormap jet;
c = colorbar;
c.Label.String = '[nm]';
view([0 90]);
axis image;
title(['Residual without tilt removed = ' num2str(nanstd(testR1(:),1)) ' nm']);


% testResult = RemoveSurface1(X,Y,Z-fitResult(32:57,:)*1e-6);
Z_fitted = RemoveSurface1(X,Y,fitResult(31:56,:)*1e3);
Z_real = RemoveSurface1(X,Y,Z*1e9);
testR2 = Z_real - Z_fitted;
subplot(414);
surf(testR2, 'EdgeColor', 'none');
% axis image;
colormap jet;
c = colorbar;
c.Label.String = '[nm]';
view([0 90]);
axis image;
title(['Residual with tilt removed = ' num2str(nanstd(testR2(:),1)) ' nm']);



    
    
