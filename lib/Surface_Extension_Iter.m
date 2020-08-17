function [TT, B, Z_residual_ca] = Surface_Extension_Iter(...
    X, Y, Z,...unextended surface error map
    brf_params,...TIF sampling interval [m/pxl]
    brf_mode,...
    du_ibf,...
    rms_thrd,...
    X_tif, Y_tif, Z_tif, ...TIF profile
    initExtMethod,...initial extention method
    ini_isFall,...
    ini_fu_range, ini_fv_range,...for gerchberg initial extention method
    ini_order_m, ini_order_n,ini_type,...  for poly initial extension method    
    iterExtMethod, ... has initial extension or not
    iter_isFall,...
    iter_fu_range, iter_fv_range,...for gerchberg initial extention method
    iter_order_m, iter_order_n,iter_type...  for poly initial extension method 
)

%% Initial extension
[X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension(...
    X,Y,Z,...
    brf_params,...
    Z_tif,...
    initExtMethod,...
    ini_isFall,...
    ini_fu_range, ini_fv_range,...
    ini_order_m, ini_order_n,ini_type);
Z_ini = Z_ext - nanmin(Z_ext(:));

%% Calculate dwell time
TT = 0;
pixel_m = median(diff(X(1,:)));
% dwell time calculation
% options = struct(...
%     'Algorithm', 'FFT',...
%     'maxIters', 10, ...
%     'PV_dif', 0.001e-9, ...[m]
%     'RMS_dif', 0.02e-9, ...[m]
%     'dwellTime_dif', 60, ...[s]
%     'isDownSampling', false, ...
%     'samplingInterval', du_ibf ... [m]
% );

options = struct(...
    'Algorithm', 'Iterative-FFT',...
    'maxIters', 50, ...
    'PV_dif', 0.01e-9, ...[m]
    'RMS_dif', 0.02e-9, ...[m]
    'dwellTime_dif', 60, ...[s]
    'isDownSampling', false, ...
    'samplingInterval', du_ibf ... [m]
);

ratio = 1;
tmin=0;
tmax=1;

[~... BRF
    , ~, ~...BRF Coordinates
    , ~, ~... full aperture results [m]
    , T_P... dwell time on the dwell grid [s]
    , ~...
    , ~, ~...dwell grid 
    , ~, ~, ~ ... dwell grid coordinates [m]
    , ~, Z_removal_dw, ~...dwell grid results [m]
    , ~, ~... clear aperture coordinates [m]
    , ~, ~, Z_residual_ca...[m]
    ] = DwellTime2D_FFT_Full_Test...
    ( X_ext, Y_ext, Z_ext... height to remove [m]
    , 0 ...
    , brf_params... BRF parameters
    , brf_mode...
    , X_tif, Y_tif, Z_tif...
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , tmin,tmax...
    , options...
    , ratio...
    , false...
    );

TT = TT + T_P;


%% Iterative refinement
Z_residual_ca_prev = Z_residual_ca;
mau_iter = 20;
iter = 1;

while(true)
  
    if(std(Z_residual_ca(:),1)<1e-9)
        [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension(...
            X,Y,Z_residual_ca,...
            brf_params,...
            Z_tif,...
            'zero', ...
            false);
    else
        [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension(...
            X,Y,Z_residual_ca,...
            brf_params,...
            Z_tif,...
            iterExtMethod,...
            iter_isFall,...
            iter_fu_range, iter_fv_range,...
            iter_order_m, iter_order_n,iter_type);
    end


    % Calculate dwell time
    du_ibf = 1e-3;     % [m] ibf dwell grid sampling interval

    % dwell time calculation
    options = struct(...
        'Algorithm', 'FFT',...
        'maxIters', 10, ...
        'PV_dif', 0.001e-9, ...[m]
        'RMS_dif', 0.02e-9, ...[m]
        'dwellTime_dif', 60, ...[s]
        'isDownSampling', false, ...
        'samplingInterval', du_ibf ... [m]
    );

    % options = struct(...
    %     'Algorithm', 'Iterative-FFT',...
    %     'maxIters', 10, ...
    %     'PV_dif', 0.001e-9, ...[m]
    %     'RMS_dif', 0.02e-9, ...[m]
    %     'dwellTime_dif', 60, ...[s]
    %     'isDownSampling', false, ...
    %     'samplingInterval', du_ibf ... [m]
    % );

    %  X_brf=0;
    %  Y_brf=0;
    %  Z_avg=0;
    ratio = 1;
    tmin=0;
    tmax=1;

    [B... BRF
        , ~, ~...BRF Coordinates
        , ~, ~... full aperture results [m]
        , T_P... dwell time on the dwell grid [s]
        , ~...
        , ~, ~...dwell grid 
        , ~, ~, ~ ... dwell grid coordinates [m]
        , ~, Z_removal_dw, ~...dwell grid results [m]
        , X_ca, Y_ca... clear aperture coordinates [m]
        , ~, ~, Z_residual_ca...[m]
        ] = DwellTime2D_FFT_Full_Test...
        ( X_ext, Y_ext, Z_ext... height to remove [m]
        , Z_removal_dw...
        , brf_params... BRF parameters
        , brf_mode...
        , X_tif, Y_tif, Z_tif...
        , ca_range... Clear aperture range [pixel] 
        , pixel_m... pixel size [m/pixel]
        , tmin,tmax...
        , options...
        , ratio...
        , false...
        );

    std_curr = std(Z_residual_ca(:),1);
    
    if(std_curr > std(Z_residual_ca_prev(:),1))
        break;
    elseif(std_curr < rms_thrd)
        TT = TT + T_P;
        break;
    elseif(iter>mau_iter)
        break;
    else
        Z_residual_ca_prev = Z_residual_ca;
        TT = TT + T_P;
        iter = iter + 1;
    end

end

Z_ext = NaN(size(X_ext));   ... mark the Z_ext to NaN
Z_ext(ca_range.v_s:ca_range.v_e, ca_range.u_s:ca_range.u_e) = Z;... fill in the valid data point
r = max(size(Z_ext)-size(Z))/2;

% obtain the are of extension
[u,v] = meshgrid(-r:r,-r:r);
[~, rr] = cart2pol(u,v);
se = rr<=r;
BW_Z = imdilate(~isnan(Z_ext), se);
id_ext = BW_Z ==0;

TT(id_ext) = 0;

Z_removal_dw = ConvFFT2(TT,B);
Z_residual_dw = Z_ini - Z_removal_dw;
Z_residual_ca = Z_residual_dw(ca_range.v_s:ca_range.v_e, ca_range.u_s:ca_range.u_e);
Z_residual_ca = RemoveSurface1(X_ca, Y_ca, Z_residual_ca);

end