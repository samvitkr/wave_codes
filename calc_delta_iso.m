clear all
close all
baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
load(fullfile(baseDir,'grid.mat'),'X','Z');
load(fullfile(baseDir,'phi_interp_2d.mat'),'uphi');
Zs = double(squeeze(Z(:,1,:)));
Xs = double(squeeze(X(:,1,:)));
pex=0.5;
Lx=2*pi/pex;
uphi = squeeze(uphi);
clear Z;
delta = uphi.*0;
[Nx, Nz]=size(Zs);

for i =1:Nx;
u=uphi(i,:);
z = Zs(i,:);
uf=double(flip(u)');
zf=double(flip(z)');

delta_i=flip(cumtrapz(zf,uf)+1);
delta_i = delta_i';
delta(i,:)=delta_i;
end


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
%% Pre-allocate your final coordinate matrices (Nz x num_contours2)
load(fullfile(baseDir,'potexact.mat')); % Loads 'phiexact'

num_contours2 = 256; % Adjust to the number of contours you want
Xphi = zeros(Nz, num_contours2);
Zphi = zeros(Nz, num_contours2);

% 1. Loop through each contour channel
for nc2 = 1:num_contours2
    
    % Enforce your starting condition at the top boundary (Nz)
    Xphi(Nz, nc2) = Xs(nc2, Nz);
    Zphi(Nz, nc2) = Zs(nc2, Nz);
    
    % Find the specific target 'phiexact' value at this top grid point
    target_val_phi = interpolateSolution(phiexact, Xphi(Nz, nc2), Zphi(Nz, nc2));
    
    % 2. Trace this specific phiexact value downwards (j = Nz-1 down to 1)
    for j = Nz-1:-1:2
        
        % Extract the X and Z coordinate vectors for the entire curved grid layer 'j'
        x_layer = Xs(:, j);
        z_layer = Zs(:, j);
        
        % Evaluate the exact solution across this curved layer
        phi_layer = interpolateSolution(phiexact, x_layer, z_layer);
        
        % Remove duplicates/non-unique values to ensure strictly monotonic points for interp1
        [phi_unique, unique_idx] = unique(phi_layer);
        x_unique = x_layer(unique_idx);
        z_unique = z_layer(unique_idx); % We must now extract unique Zs as well
        
        % 3. Interpolate to find BOTH the exact X and Z position where phi == target_val_phi
        if length(phi_unique) > 1 && ...
           target_val_phi >= min(phi_unique) && ...
           target_val_phi <= max(phi_unique)
            
            % Interpolate along the curve to find the exact (X, Z) coordinate
            Xphi(j, nc2) = interp1(phi_unique, x_unique, target_val_phi, 'linear');
            Zphi(j, nc2) = interp1(phi_unique, z_unique, target_val_phi, 'linear');
        else
            % If the target value leaves the boundary of this layer, pad with NaN
            Xphi(j, nc2) = NaN;
            Zphi(j, nc2) = NaN;
        end
    end
end
%%
X_phi = Xphi(2:end,:)';
Z_phi = Zphi(2:end,:)';
save(fullfile(baseDir,"slines.mat"),'Zsl','X_phi','Z_phi','-append')

%% Calculate Intersections using custom 'intersections.m'

% Pre-allocate the orthogonal grid matrices with NaNs
% Dim 1 (i): constant phi contour
% Dim 2 (j): constant stream function contour
Xgphi = NaN(num_contours2, num_contours);
Zgphi = NaN(num_contours2, num_contours);

% 1. Loop over each constant phi contour (index 'i')
for i = 1:num_contours2
    
    % Extract the current phi contour (column 'i')
    x_phi = Xphi(2:end, i);
    z_phi = Zphi(2:end, i);
    
    % Strip NaNs to optimize the bounding box checks inside intersections.m
    valid_phi = ~isnan(x_phi) & ~isnan(z_phi);
    x_phi = x_phi(valid_phi);
    z_phi = z_phi(valid_phi);
    
    % 2. Loop over each constant stream function contour (index 'j')
    for j = 1:num_contours
        
        % Extract the current stream function contour (column 'j')
        x_psi = Xsl(:, j);
        z_psi = Zsl(:, j);
        
        % Strip NaNs
        valid_psi = ~isnan(x_psi) & ~isnan(z_psi);
        x_psi = x_psi(valid_psi);
        z_psi = z_psi(valid_psi);
        
        % Proceed only if both curves have at least 2 points to form a segment
        if length(x_phi) > 1 && length(x_psi) > 1
            
            % Call the uploaded intersections function
            [x_int, z_int] = intersections(x_phi, z_phi, x_psi, z_psi);
            
            % If one or more intersections exist, store the primary one
            if ~isempty(x_int)
                Xgphi(i, j) = x_int(1);
                Zgphi(i, j) = z_int(1);
            end
            
        end
    end
end

%%
Xgphi=Xgphi(:,2:end);
Zgphi=Zgphi(:,2:end);
save(fullfile(baseDir,'uniform_phi_grid.mat'),'Xgphi','Zgphi','-v7.3')
% Optional: Visualize the final generated grid
%figure;
%hold on;
%% Plot the phi contours in red
%plot(Xphi, Zphi, 'r-', 'LineWidth', 0.5);
%% Plot the stream function contours in blue
%plot(Xsl, Zsl, 'b-', 'LineWidth', 0.5);
%% Plot the intersection points in black
%plot(Xgphi(:), Zgphi(:), 'k.', 'MarkerSize', 10);
%hold off;
%title('Orthogonal Grid Intersections');
%xlabel('X');
%ylabel('Z');
%axis equal;
