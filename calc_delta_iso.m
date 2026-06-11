clear all
close all
baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
load(fullfile(baseDir,'grid.mat'),'X','Z');
load(fullfile(baseDir,'phi_interp_2d.mat'),'uphi');
Zs = squeeze(Z(:,1,:));
Xs = squeeze(X(:,1,:));
pex=0.5;
Lx=2*pi/pex;
uphi = squeeze(uphi);
clear Z;
delta = uphi.*0;
[Nx, Nz]=size(Zs);

i=1;
for i =1:Nx;
u=uphi(1,:);
z = Zs(1,:);
uf=double(flip(u)');
zf=double(flip(z)');

delta_i=flip(cumtrapz(zf,uf)+1);
delta_i = delta_i';
delta(i,:)=delta_i;
end

% 1. Initialize your output matrices (Nx rows by 120 columns)
num_contours = 128;
Xsl = zeros(Nx, num_contours);
Zsl = zeros(Nx, num_contours);

% 2. Loop through each contour channel
for nc = 1:num_contours
    
    % Enforce your starting condition explicitly in the first row
    Xsl(1, nc) = Xs(1, nc);
    Zsl(1, nc) = Zs(1, nc);
    
    % Find the specific target 'delta' value at this starting grid point
    target_val = delta(1, nc);
    
    % 3. Trace this specific delta value across the remaining rows (i = 2 to Nx)
    for i = 2:Nx
        % Assign the standard X coordinate for this row/contour profile
        Xsl(i, nc) = Xs(i, nc);
        
        % Extract the delta and Z profiles for the current row 'i'
        delta_row = delta(i, :);
        z_row     = Zs(i, :);
        
        % Remove duplicates/non-unique values to ensure clean interpolation
        [delta_unique, unique_idx] = unique(delta_row);
        z_unique = z_row(unique_idx);
        
        % 4. Interpolate to find the exact Z position where delta == target_val
        if length(delta_unique) > 1 && ...
           target_val >= min(delta_unique) && ...
           target_val <= max(delta_unique)
       
            Zsl(i, nc) = interp1(delta_unique, z_unique, target_val, 'linear');
        else
            % If the target value leaves the boundary of this row, pad with NaN
            Zsl(i, nc) = NaN; 
        end
    end
end
Zdelta = Zsl;
%%
x_uniform_phi = zeros(Nx, num_contours);
phi_target = Lx * (0:Nx-1)' / Nx;
 load(fullfile(baseDir,'potexact.mat'));

for idx = 1:num_contours

    z_line = double(Zsl(:, idx));
    x_1D = double(Xs(:,idx));
    phi_current = interpolateSolution(phiexact, x_1D, z_line);
    phi_current = Lx * (phi_current - phi_current(1)) / (phi_current(end) - phi_current(1));
    [phi_clean, unique_idx] = unique(phi_current);
    x_clean = x_1D(unique_idx);
    x_uniform_phi(:, idx) = interp1(phi_clean, x_clean, phi_target, 'spline', 'extrap');
    %slq(:,idx)=cos(k0 * x_uniform_phi(:,idx)) * asl(idx) + msl(idx);
end
xq=x_uniform_phi;

%sl = fullfile(baseDir,'slines.mat');
%save(sl,'Zdelta','delta','-append')
