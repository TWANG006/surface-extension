function Z_fall = Surface_Extension_Fall(...
    X_ext, Y_ext, Z_ext,...extended surface
    ca_range,...ca range in pixels
    brf_params,...brf parameters
    Z_tif...tif profile
    )
% Function
%   Z_fall = Surface_Extension_Fall(X_ext, Y_ext, Z_ext, ca_range,...
%                                   tif_mpp, Z_tif)
% Purpose
%   Apply the fall profile to the extended part of the surface

%% Obtain parameters
r = max(size(Z_tif)-1) * brf_params.lat_res_brf * 0.5; ... radius of the TIF

Z_fall = NaN(size(Z_ext));   ... initial extension matrix with NaNs
Z_fall(ca_range.v_s:ca_range.v_e, ca_range.u_s:ca_range.u_e) = 1;    ... fill 1 in the valid data point

%% Finding edge points
id_edg = Surface_Extension_EdgeExtraction(Z_fall);
u_edg = X_ext(id_edg);
v_edg = Y_ext(id_edg);

% figure,imagesc(id_edg);

%% Obtain the filled & original ids
id_fil = isnan(Z_fall); ... filled data ids
u_fil = X_ext(id_fil);  ... x coordinates of filled data
v_fil = Y_ext(id_fil);  ... y coordinates of filled data

%% Calculate fall profiles
fun = @(x, A, sigma) A*exp(-(x).^2/(2*sigma.^2));
B = 1/integral(@(x)fun(x,brf_params.A, brf_params.sigma_xy(1)), -(r), r);

fall_profiles = zeros*u_fil;
for k = 1:length(u_fil)
    % min distances from filled points to edge points
%     fall_profiles(k) = min(sqrt((u_fil(k) - u_edg).^2+(v_fil(k) - v_edg).^2));
    
    % calculate the fall profile
    fall_profiles(k) = B*integral(@(x)fun(x,brf_params.A, brf_params.sigma_xy(1)), -(r-min(sqrt((u_fil(k) - u_edg).^2+(v_fil(k) - v_edg).^2))), r);
    
end

Z_fall(id_fil) = fall_profiles;
Z_fall = Z_ext.*Z_fall;

% fsfig('');
% subplot(311);
% surf(X_ext, Y_ext, Z_ext, 'EdgeColor', 'none');
% % rectangle('Position', [xs, ys, (je-js+1), (ie-is+1)]);
% axis image;
% view([0 90]);
% 
% subplot(312);
% surf(X_ext, Y_ext, Z_fall, 'EdgeColor', 'none');
% % rectangle('Position', [xs, ys, (je-js+1), (ie-is+1)]);
% axis image;
% view([0 90]);
% 
% 
% 
% subplot(313);
% surf(X_ext, Y_ext, Z_fall, 'EdgeColor', 'none');
% % rectangle('Position', [xs, ys, (je-js+1), (ie-is+1)]);
% axis image;
% view([0 90]);



end