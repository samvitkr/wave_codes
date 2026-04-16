%% Load grid
data = load('grid.mat');    % expects X, Y, Zw
X  = data.X;
Y  = data.Y;
Zw = data.Zw;

[Nx,Ny,Nz] = size(Zw);

% Take one representative Y-slice (solution is uniform in y)
iy = 1;
x  = squeeze(X(:,iy,1));      % Nx x 1
z  = squeeze(Zw(:,iy,:));     % Nx x Nz  (each row: varying z along column)

% Sanity checks
if any(diff(x) <= 0)
    error('x must be strictly increasing.');
end
if size(z,1) ~= Nx || size(z,2) ~= Nz
    error('Unexpected z size; check Zw layout.');
end

%% Indices and spacings
N  = Nx * Nz;
A  = spalloc(N, N, 5*N);
b  = zeros(N,1);

idx = @(i,k) (k-1)*Nx + i;    % map (i,k) -> 1D index

dx = diff(x);                 % length Nx-1

% z spacings
dz_forward  = zeros(Nx,Nz);
dz_backward = zeros(Nx,Nz);
for i = 1:Nx
    zi = z(i,:);
    dz_forward(i,1:Nz-1) = zi(2:end) - zi(1:end-1);
    dz_forward(i,Nz)     = dz_forward(i,Nz-1);
    dz_backward(i,2:Nz)  = zi(2:end) - zi(1:end-1);
    dz_backward(i,1)     = dz_backward(i,2);
end

%% Assemble A, b
for k = 1:Nz
    for i = 1:Nx
        p = idx(i,k);

        % Dirichlet at inflow/outflow in x
        if i == 1
            A(p,p) = 1.0;
            b(p)   = -0.5;
            continue;
        elseif i == Nx
            A(p,p) = 1.0;
            b(p)   = 0.5;
            continue;
        end

        % Neumann at bottom and top in z (zero normal gradient)
        if k == 1
            % bottom: phi(i,2) - phi(i,1) = 0
            A(p,p)        = -1.0;
            A(p,idx(i,2)) =  1.0;
            b(p)          =  0.0;
            continue;
        elseif k == Nz
            % top: phi(i,Nz) - phi(i,Nz-1) = 0
            A(p,p)           =  1.0;
            A(p,idx(i,Nz-1)) = -1.0;
            b(p)             =  0.0;
            continue;
        end

        % Interior: non-uniform Laplacian in x,z
        % x second derivative (1D non-uniform)
        dxm = dx(i-1);          % x(i)   - x(i-1)
        dxp = dx(i);            % x(i+1) - x(i)
        axm =  2 / (dxm*(dxm+dxp));
        axp =  2 / (dxp*(dxm+dxp));
        ax0 = -axm - axp;

        % z second derivative (non-uniform, per i)
        dzm = dz_backward(i,k); % z(i,k) - z(i,k-1)
        dzp = dz_forward(i,k);  % z(i,k+1) - z(i,k)
        azm =  2 / (dzm*(dzm+dzp));
        azp =  2 / (dzp*(dzm+dzp));
        az0 = -azm - azp;

        A(p,idx(i-1,k)) = axm;
        A(p,idx(i+1,k)) = axp;
        A(p,idx(i,k-1)) = azm;
        A(p,idx(i,k+1)) = azp;
        A(p,p)          = ax0 + az0;
        b(p)            = 0.0;
    end
end

%% Solve
% Optional: quick diagnostic
% figure; spy(A); title('Sparsity of A');

phi_vec = A \ b;

% Reshape to x-z plane
phi_xz = reshape(phi_vec,[Nx,Nz]);   % 128 x 128

%% Extend uniformly in y
phi = repmat(phi_xz, [1, Ny, 1]);    % 128 x 16 x 128

%% Gradients in physical space
dphidx = zeros(Nx,Nz);
dphidz = zeros(Nx,Nz);

% d/dx, non-uniform
for k = 1:Nz
    % interior
    for i = 2:Nx-1
        dxm = x(i)   - x(i-1);
        dxp = x(i+1) - x(i);
        dphidx(i,k) = (-(dxp)/(dxm*(dxm+dxp))*phi_xz(i-1,k) ...
                       + (dxp-dxm)/(dxm*dxp)*phi_xz(i,k) ...
                       + (dxm)/(dxp*(dxm+dxp))*phi_xz(i+1,k));
    end
    % boundaries
    dphidx(1,k)  = (phi_xz(2,k)    - phi_xz(1,k))   / (x(2)-x(1));
    dphidx(Nx,k) = (phi_xz(Nx,k)   - phi_xz(Nx-1,k))/(x(Nx)-x(Nx-1));
end

% d/dz, non-uniform per i
for i = 1:Nx
    zi = z(i,:);
    % interior
    for k = 2:Nz-1
        dzm = zi(k)   - zi(k-1);
        dzp = zi(k+1) - zi(k);
        dphidz(i,k) = (-(dzp)/(dzm*(dzm+dzp))*phi_xz(i,k-1) ...
                       + (dzp-dzm)/(dzm*dzp)*phi_xz(i,k) ...
                       + (dzm)/(dzp*(dzm+dzp))*phi_xz(i,k+1));
    end
    % boundaries (consistent with zero-normal BC)
    dphidz(i,1)  = (phi_xz(i,2)    - phi_xz(i,1))   / (zi(2)-zi(1));
    dphidz(i,Nz) = (phi_xz(i,Nz)   - phi_xz(i,Nz-1))/(zi(Nz)-zi(Nz-1));
end

dphidx_3d = repmat(dphidx, [1, Ny, 1]);
dphidy_3d = zeros(Nx, Ny, Nz);
dphidz_3d = repmat(dphidz, [1, Ny, 1]);

%% Plot solution (x-z plane at chosen y-slice)
figure;
surf(squeeze(X(:,iy,:)), squeeze(Zw(:,iy,:)), squeeze(phi(:,iy,:)));
shading interp; colorbar;
xlabel('x'); ylabel('z'); zlabel('\phi');
title(sprintf('Laplace solution at y-index %d', iy));
view(2);    % comment this out if you prefer 3D view

