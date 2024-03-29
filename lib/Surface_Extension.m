function [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension(...
    X, Y, Z,...unextended surface error map
    brf_params,...TIF sampling interval [m/pxl]
    Z_tif,  ...TIF profile
    method, ...extension method
    isFall, ...using fall profile or not
    fu_range, fv_range, ...frequency domains for Gerchberg Pouplis algorithm
    order_m, order_n,... orders for polynomial fitting algorithm
    type...polynomial type for polynomial fitting, 'Chebyshev' or 'Legendre'
    )
%% Default parameters
if nargin == 5
    method = 'smooth';
    isFall = false;
end
if nargin == 6
    isFall = false;
end

%% Different extension algorithms
if strcmp(method, 'zero')
    [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension_Zeros(X,Y,Z,brf_params.lat_res_brf,Z_tif);
elseif strcmp(method, '8nn')
    [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension_8NN(X,Y,Z,brf_params.lat_res_brf,Z_tif);
elseif strcmp(method, 'smooth')
    [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension_Smooth(X,Y,Z,brf_params.lat_res_brf,Z_tif);
elseif strcmp(method, 'gauss')
    [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension_Gauss(X,Y,Z,brf_params,Z_tif);
elseif strcmp(method, 'gerchberg')
    if nargin ~= 9
        error('Not enough parameters for Gerchberg algorithm:fu_range and fv_range should be fed.');
    else
        [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension_GP(X,Y,Z,brf_params.lat_res_brf,Z_tif, fu_range, fv_range);
    end
elseif strcmp(method, 'poly')
    if nargin ~=12
        error('Not enough parameters for polynomial fitting algorithm: order_m, order_n, and type should be fed.');
    else
        [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension_Polyfit(X,Y,Z,brf_params.lat_res_brf,Z_tif,order_m,order_n,type);
    end
else
    error('Invalid algorithm selected.');
end

%% Fall or not
if isFall
    if(~strcmp(method, 'zeros') && ~strcmp(method, 'gauss') && ~strcmp(method, 'poly'))
        Z_ext = Surface_Extension_Fall(X_ext, Y_ext, Z_ext, ca_range, brf_params, Z_tif);
    else
        warning(['Fall profile is automatically disabled for ' method ' algorithm']);
    end
end

end