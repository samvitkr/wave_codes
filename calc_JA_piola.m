

clear;
close all;

% ==========================================
% 1. PARAMETERS & INITIALIZATION
% ==========================================
baseDir = '/scratch.global/kuma0458/c8ak1_re180/run';
c       = 8;
ak 	= 0.1;
ret     = 180;
nu      = 1 / ret;
Nx      = 256;
Ny      = 192;
Nz      = 128;

tstart  = 4300000000;
step    =    1250000;
tend    = 4770000000;

% ==========================================
% 2. LOAD STATIC GRID DEFINITIONS
% ==========================================
fname_gridh5 = fullfile(baseDir, 'grid.h5');
fprintf('Loading computational grid definitions...\n');
zz  = h5read(fname_gridh5, '/zz');
zw  = h5read(fname_gridh5, '/zw');
pex = h5read(fname_gridh5, '/pex');
pey = h5read(fname_gridh5, '/pey');

Lx = 2 * pi / pex;
Ly = 2 * pi / pey;
H  = 1; % Mean height

% Wavenumbers for Fourier Derivatives
wave_n = 12;
k0     = wave_n * pex;
a_wave = ak / k0;

kx    = pex * [0:Nx/2-1, -Nx/2:-1]';
ky    = pey * [0:Ny/2-1, -Ny/2:-1]';
kx_3D = reshape(kx, [Nx, 1, 1]);
ky_3D = reshape(ky, [1, Ny, 1]);

% 1D coordinate vectors
x_1d = (0:Nx-1)' * Lx / Nx;

% Pre-allocate Time Series Arrays
times     = [];
Ia_series = [];
T_dt_a    = [];
T_Vn      = [];
T_R       = [];

% Finite Difference Stencils for Z-derivatives
dz = diff(zz);
nz = Nz;
idx_i  = 1:(nz-2);
a_st   = reshape(1./dz(idx_i+1) - 1./(dz(idx_i) + dz(idx_i+1)), 1, 1, []);
b_st   = reshape(1./dz(idx_i)   - 1./(dz(idx_i+1)), 1, 1, []);
c_st   = reshape(1./(dz(idx_i) + dz(idx_i+1))  - 1./dz(idx_i), 1, 1, []);

% ==========================================
% 3. MAIN INTEGRATION LOOP
% ==========================================
for tstep = tstart:step:tend
    % Load Velocity Fields (Lab Frame)
    fname = fullfile(baseDir, sprintf('Sol%014d.h5', tstep));
    fprintf('Processing %s\n', fname);

    u    = h5read(fname, '/u');
    v    = h5read(fname, '/v');
    w    = h5read(fname, '/w');
    time = h5read(fname, '/time');
    times = [times; time];

    % Interpolate w to cell centers (computational space z-direction)
    wc = permute(w, [3 1 2]);
    wc = interp1(zw(1:end-1), wc(1:end-1, :, :), zz);
    wc = single(permute(wc, [2 3 1]));

    % Load Moving Grid Fields at current time step[cite: 3]
    fng = fullfile(baseDir, sprintf('grid%014d.mat', tstep));
    load(fng, 'Z', 'dZetadz', 'dZetadx');

    % Ensure correct metric dimensions[cite: 3]
    dZetadz_3D = reshape(dZetadz, [Nx, 1, 1]);
    dZetadx_3D = reshape(dZetadx, [Nx, 1, Nz]);
    Jac_1D     = squeeze(1 ./ dZetadz(:,1)); % Jacobian J = 1/dZetadz

    % --- A. PIOLA CARRIER & EULERIAN TIME DERIVATIVES ---
    % Since a_x = H/d and a_z = H(H-z)\eta_x / d^2[cite: 9]
    % It corresponds exactly to the moving metrics[cite: 3]
    a_x_3D = single(H .* dZetadz_3D);
    a_z_3D = single(-H .* dZetadx_3D);

    % Analytical time derivatives for precision (evaluated at lab-frame x)[cite: 9]
    phase   = k0 * (x_1d - c * time);
    eta     = a_wave * cos(phase);
    eta_x   = -a_wave * k0 * sin(phase);
    eta_xx  = -a_wave * k0^2 * cos(phase);

    eta_t   = -c * eta_x;
    eta_xt  = -c * eta_xx;
    d_gap   = H - eta;

    % \partial_t a evaluated at instantaneous Z physical height[cite: 9]
    dt_a_x = H .* eta_t ./ (d_gap.^2);
    dt_a_z = H .* (H - squeeze(Z(:,1,:))) .* ( (eta_xt ./ (d_gap.^2)) + (2 .* eta_x .* eta_t ./ (d_gap.^3)) );

    dt_a_x_3D = reshape(dt_a_x, [Nx, 1, 1]);
    dt_a_z_3D = reshape(dt_a_z, [Nx, 1, Nz]);

    % --- B. COMPUTE SPATIAL GRADIENTS (LAB FRAME) ---
    % Vertical (zeta) derivatives
    dudzeta = zeros(size(u), 'single');
    dvdzeta = zeros(size(v), 'single');
    dwdzeta = zeros(size(wc), 'single');

    dudzeta(:,:,1) = (u(:,:,2) - u(:,:,1)) ./ dz(1);
    dvdzeta(:,:,1) = (v(:,:,2) - v(:,:,1)) ./ dz(1);
    dwdzeta(:,:,1) = (wc(:,:,2) - wc(:,:,1)) ./ dz(1);

    dudzeta(:,:,end) = (u(:,:,end) - u(:,:,end-1)) ./ dz(end);
    dvdzeta(:,:,end) = (v(:,:,end) - v(:,:,end-1)) ./ dz(end);
    dwdzeta(:,:,end) = (wc(:,:,end)- wc(:,:,end-1)) ./ dz(end);

    dudzeta(:,:,2:nz-1) = a_st.*u(:,:,3:nz)  + b_st.*u(:,:,2:nz-1)  + c_st.*u(:,:,1:nz-2);
    dvdzeta(:,:,2:nz-1) = a_st.*v(:,:,3:nz)  + b_st.*v(:,:,2:nz-1)  + c_st.*v(:,:,1:nz-2);
    dwdzeta(:,:,2:nz-1) = a_st.*wc(:,:,3:nz) + b_st.*wc(:,:,2:nz-1) + c_st.*wc(:,:,1:nz-2);

    dudz = single(dZetadz_3D .* dudzeta);
    dvdz = single(dZetadz_3D .* dvdzeta);
    dwdz = single(dZetadz_3D .* dwdzeta);

    % Horizontal Fourier Derivatives (Direct Lab Frame computation)
    dudy = single(ifft(fft(u, [], 2)  .* (1i .* ky_3D), [], 2, 'symmetric'));
    dvdy = single(ifft(fft(v, [], 2)  .* (1i .* ky_3D), [], 2, 'symmetric'));
    dwdy = single(ifft(fft(wc, [], 2) .* (1i .* ky_3D), [], 2, 'symmetric'));

    dudx = single(ifft(fft(u, [], 1)  .* (1i .* kx_3D), [], 1, 'symmetric') + dZetadx_3D .* dudzeta);
    dvdx = single(ifft(fft(v, [], 1)  .* (1i .* kx_3D), [], 1, 'symmetric') + dZetadx_3D .* dvdzeta);
    dwdx = single(ifft(fft(wc, [], 1) .* (1i .* kx_3D), [], 1, 'symmetric') + dZetadx_3D .* dwdzeta);

    % Viscous Laplacians
    Jinv = 1 ./ dZetadz_3D;

    F1_u = Jinv .* dudx;
    dF1_u_dxi = ifft(fft(F1_u, [], 1) .* (1i .* kx_3D), [], 1, 'symmetric');
    F3_u = Jinv .* (dZetadx_3D .* dudx + (dZetadz_3D.^2) .* dudzeta);
    dF3_u_dzeta = zeros(size(u), 'single');
    dF3_u_dzeta(:,:,1) = (F3_u(:,:,2) - F3_u(:,:,1)) ./ dz(1);
    dF3_u_dzeta(:,:,2:nz-1) = a_st.*F3_u(:,:,3:nz) + b_st.*F3_u(:,:,2:nz-1) + c_st.*F3_u(:,:,1:nz-2);
    dF3_u_dzeta(:,:,end) = (F3_u(:,:,end) - F3_u(:,:,end-1)) ./ dz(end);
    viscu = single(nu .* ((dF1_u_dxi + dF3_u_dzeta) ./ Jinv + ifft(fft(u,[],2).*(-(ky_3D).^2),[],2,'symmetric')));

    F1_w = Jinv .* dwdx;
    dF1_w_dxi = ifft(fft(F1_w, [], 1) .* (1i .* kx_3D), [], 1, 'symmetric');
    F3_w = Jinv .* (dZetadx_3D .* dwdx + (dZetadz_3D.^2) .* dwdzeta);
    dF3_w_dzeta = zeros(size(wc), 'single');
    dF3_w_dzeta(:,:,1) = (F3_w(:,:,2) - F3_w(:,:,1)) ./ dz(1);
    dF3_w_dzeta(:,:,2:nz-1) = a_st.*F3_w(:,:,3:nz) + b_st.*F3_w(:,:,2:nz-1) + c_st.*F3_w(:,:,1:nz-2);
    dF3_w_dzeta(:,:,end) = (F3_w(:,:,end) - F3_w(:,:,end-1)) ./ dz(end);
    viscw = single(nu .* ((dF1_w_dxi + dF3_w_dzeta) ./ Jinv + ifft(fft(wc,[],2).*(-(ky_3D).^2),[],2,'symmetric')));

    % --- C. VORTICITY & FLUX TENSOR (R) ---
    ox = dwdy - dvdz;
    oy = dudz - dwdx;
    oz = dvdx - dudy;

    % R = u x w - nu * curl(w)[cite: 9]
    R_x = v .* oz - wc .* oy + viscu;
    R_z = u .* oy - v .* ox + viscw;

    % --- D. VOLUME & SURFACE INTEGRALS FOR RELATION 8 ---

    % Integral 1: Carrier Momentum (I_a)[cite: 9]
    Ia_field = a_x_3D .* u + a_z_3D .* wc;
    Ia_mean  = squeeze(mean(Ia_field, 2)); % Average over span
    Ia       = trapz(x_1d, Jac_1D .* trapz(zz, Ia_mean, 2));
    Ia_series = [Ia_series; Ia];

    % Integral 2: Eulerian Change (u * dt_a)[cite: 9]
    udta_field = u .* dt_a_x_3D + wc .* dt_a_z_3D;
    udta_mean  = squeeze(mean(udta_field, 2));
    T_dt_a     = [T_dt_a; trapz(x_1d, Jac_1D .* trapz(zz, udta_mean, 2))];

    % Integral 3: Vorticity Flux (a * R)[cite: 9]
    aR_field = a_x_3D .* R_x + a_z_3D .* R_z;
    aR_mean  = squeeze(mean(aR_field, 2));
    T_R      = [T_R; trapz(x_1d, Jac_1D .* trapz(zz, aR_mean, 2))];

    % Integral 4: Moving Surface Flux (Vn a * u)[cite: 9]
    % V_n dS = -\eta_t dx dy[cite: 9]
    u_wall    = squeeze(mean(u(:,:,1), 2));
    w_wall    = squeeze(mean(wc(:,:,1), 2));
    ax_wall   = squeeze(a_x_3D(:,:,1));
    az_wall   = squeeze(a_z_3D(:,:,1));
	Vn = -eta_t./(sqrt(1 + eta_x.^2));
    surface_flux = Vn.* (ax_wall .* u_wall + az_wall .* w_wall);
    T_Vn         = [T_Vn; trapz(x_1d, surface_flux)];
end

% ==========================================
% 4. COMPUTE TIME DERIVATIVE & VERIFY BALANCE
% ==========================================
dIa_dt = zeros(size(Ia_series));
if length(times) > 2
    dIa_dt(2:end-1) = (Ia_series(3:end) - Ia_series(1:end-2)) ./ (times(3:end) - times(1:end-2));
    dIa_dt(1)       = (Ia_series(2) - Ia_series(1)) / (times(2) - times(1));
    dIa_dt(end)     = (Ia_series(end) - Ia_series(end-1)) / (times(end) - times(end-1));
end

% Construct Relation 8 Right Hand Side: dIa/dt - Int(u*dt_a) - Int(Vn*a*u) - Int(a*R)[cite: 9]
RHS_total = dIa_dt - T_dt_a - T_Vn - T_R;

% Save results
save(fullfile(baseDir, 'Piola_Relation8_LabFrame.mat'), 'times', 'dIa_dt', 'T_dt_a', 'T_Vn', 'T_R', 'RHS_total');

% % Plot Diagnostic Balance
%figure('Name', 'Piola JA Relation 8 Balance (Lab Frame)');
%plot(times, dIa_dt, 'LineWidth', 1.5, 'DisplayName', 'dI_a/dt'); hold on;
%plot(times, T_dt_a, 'LineWidth', 1.5, 'DisplayName', '\int u \cdot \partial_t a dV');
%plot(times, T_Vn, 'LineWidth', 1.5, 'DisplayName', '\int V_n a \cdot u dS');
%plot(times, T_R, 'LineWidth', 1.5, 'DisplayName', '\int a \cdot R dV');
%plot(times, RHS_total, 'k--', 'LineWidth', 2, 'DisplayName', 'Total RHS (K \Delta B)');
%xlabel('Time');
%ylabel('Virtual Work / Flux');
%title('Piola Carrier Moving-Domain Diagnostic Balance');
%legend('Location', 'best');
%grid on;

fprintf('Relation 8 tested successfully in the laboratory frame. Diagnostic saved.\n');



