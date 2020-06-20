function [X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension_8NN(...
    X, Y, Z,...unextended surface error map
    tif_mpp,...TIF sampling interval [m/pxl]
    Z_tif...TIF profile
    )
% Function
%   [X_ext, Y_ext, Z_ext] = Surface_Extenstion_8NN(X, Y, Z, X_tif, Y_tif, Z_tif)
% Purpose
%   Extend the surface error map using 8 nearest neighbors

%% 0. Obtain required parameters
% Sampling intervals
surf_mpp = median(diff(X(1,:)));    ... surface sampling interval [m/pxl]

m = size(Z,1);  ... CA height [pixel]
n = size(Z,2);  ... CA width [pixel]

m_ext = round(tif_mpp*(size(Z_tif, 1))*0.5/surf_mpp);   ... extension size in y [pixel]
n_ext = round(tif_mpp*(size(Z_tif, 2))*0.5/surf_mpp);   ... extension size in x [pixel] 

ca_range.y_s = m_ext + 1;   ca_range.y_e = ca_range.y_s + m - 1;   ... y start & end ids of CA in FA [pixel]
ca_range.x_s = n_ext + 1;   ca_range.x_e = ca_range.x_s + n - 1;   ... x start & end ids of CA in FA [pixel]


%% 1. Initial extension matrices
[X_ext, Y_ext] = meshgrid(-n_ext:n-1+n_ext, -m_ext:m-1+m_ext);  ...extension grid
X_ext = X_ext * surf_mpp + X(1,1);  ... adjust X grid add X(1,1)
Y_ext = Y_ext * surf_mpp + Y(1,1);  ... adjust Y grid add Y(1,1)

Z_ext = NaN(size(X_ext));   ... mark the Z_ext to NaN
Z_ext(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e) = Z;... fill in the valid data points
BW_ini = ~isnan(Z_ext);   ... obtain the black&white map
BW_prev = BW_ini;
% BW_all=0;
h = size(Z_ext, 1);
w = size(Z_ext, 2);

%% 2. Filling the invalid points
r = 1;  ... extension radius (1 ~ max(m_ext, n_ext))
while(r <= max(m_ext, n_ext))

[u,v] = meshgrid(-r:r,-r:r);
[~, rr] = cart2pol(u,v);
se = rr<=r;

BW_curr = imdilate(BW_ini, se);
BW_fill = BW_curr - BW_prev;
[idy, idx] = find(BW_fill==1);

while(~isempty(idy))
    % 8-neighbor averaging
    for k = 1:length(idy)        
        count = 0;
        nn_sum = 0;

        for i = -1:1
            for j = -1:1
                if (~(i==0 && j==0))                    
                    idi = idy(k)+i;    ... neighbor y id
                    idj = idx(k)+j;    ... neighbor x id
                    
                    if (0<idi && idi<=h && 0<idj && idj<=w && ~isnan(Z_ext(idi, idj)))
                        count = count+1;
                        nn_sum = nn_sum + Z_ext(idi, idj);
                    end
                end
            end
        end

        if (count >=3)
            Z_ext(idy(k), idx(k)) = nn_sum/count;
            BW_fill(idy(k), idx(k)) = 0;
        end

    end
    [idy, idx] = find(BW_fill==1);
end

BW_prev = BW_curr;
r = r + 1;

end

Z_ext(isnan(Z_ext)) = 0;

end
    