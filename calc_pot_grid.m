clear; 
close all;
baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
fn_grid = fullfile(baseDir, 'grid.mat');
fn_pot  = fullfile(baseDir, 'potexact.mat');
fn_out  = fullfile(baseDir, 'phi_interp_2d.mat'); % Saving as a 2D specific file
fname = fullfile(baseDir,'grid.h5');
fprintf('Loading grid and exact potential...\n');
% Load X and Z (we only need the coordinates)
load(fn_grid); 
% Load the PDE stationary results object
load(fn_pot, 'phiexact'); 

[Nx, ~, Nz] = size(X);

% Squeeze out the first y-index to get a pure [Nx, Nz] matrix
X_2d = squeeze(X(:, 1, :));
Z_2d = squeeze(Z(:, 1, :));

% Flatten arrays for the PDE interpolator
xq = double(X_2d(:))';
zq = double(Z_2d(:))';

fprintf('Interpolating 2D slice using PDE toolbox...\n');
phi_2d_flat = interpolateSolution(phiexact, xq, zq);

nan_idx = isnan(phi_2d_flat);
if any(nan_idx)
    fprintf('Found %d NaN points near boundaries. Fixing with nearest-neighbor extrapolation...\n', sum(nan_idx));
    % Extract unstructured nodes and nodal solutions
    nodes = phiexact.Mesh.Nodes;
    sol = phiexact.NodalSolution;
    % Create a scattered interpolant. 
    % 'linear' interpolation inside, 'nearest' extrapolation outside prevents NaNs
    F = scatteredInterpolant(nodes(1,:)', nodes(2,:)', sol, 'linear', 'nearest');
    % Evaluate F only at the NaN locations
    phi_2d_flat(nan_idx) = F(xq(nan_idx)', zq(nan_idx)');
end

% Reshape the cleaned 1D array back into the 2D [Nx, Nz] matrix
Phi_2D = single(reshape(phi_2d_flat, [Nx, Nz])); 

fprintf('Saving 2D interpolated potential to %s...\n', fn_out);
fprintf('Interpolation complete.\n');

%%
zz   = h5read(fname, '/zz');
zw   = h5read(fname,'/zw');
pex  = h5read(fname, '/pex');

u=reshape(Phi_2D,Nx, 1, Nz)-X(:,1,:);
kx=pex*[0:Nx/2-1,-Nx/2:-1]';

dudzeta=zeros(size(u));
dz=diff(zz);
dudzeta(:,:,1)=(u(:,:,2)-u(:,:,1))./dz(1);
nz=Nz;
i      = 1:(nz-2);
dz_i   = dz(i);
dz_ip1 = dz(i+1);
sumdz  = dz_i + dz_ip1;
a  = 1./dz_ip1 - 1./sumdz;
b  = 1./dz_i   - 1./dz_ip1;
c  = 1./sumdz  - 1./dz_i;
% reshape for implicit expansion along 3rd dim
a  = reshape(a, 1, 1, []);
b  = reshape(b, 1, 1, []);
c  = reshape(c, 1, 1, []);
% Apply to u, v
dudzeta(:,:,2:nz-1) =a.*u(:,:,3:nz)+b.*u(:,:,2:nz-1)+c.* u(:,:,1:nz-2);
dudz = single(dZetadz.*dudzeta);
fxu =fft(u,[],1);
dudxi=ifft( fxu.*(1i.*kx),[],1,'symmetric');
dudx=single(dudxi+dZetadx.*(dudzeta));

uphi=dudx+1.0;
wphi=dudz;

dzw =diff(zw)';
Jacobian=1./dZetadz;
uslice = reshape(uphi,Nx,Nz);
uslice(isnan(uslice))=0;
udz=sum(uslice(:,2:end).*dzw,2);
Jx=udz.*Jacobian;
J=mean(Jx);
uphi=uphi./J;
Phi_2D=Phi_2D./J;
wphi=wphi./J;

save(fn_out, 'J','Phi_2D','uphi','wphi', '-v7.3');
