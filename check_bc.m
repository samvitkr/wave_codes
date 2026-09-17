

clear; close all;

% Define base parameters (adjust these matching your simulation run)
baseDir = '/scratch.global/kuma0458/c8ak1_re180/run';
tstep = 4770000000; % Replace with specific time step

c = 8;
ak = 0.1;

Nx = 256;
Ny = 192;
wave_n=12;
% Construct filenames
gridName = fullfile(baseDir, 'grid.h5');
% Formats the filename with 15 zero-padded digits as requested
fname = fullfile(baseDir, sprintf('Sol%014d.h5', tstep));

% Read grid parameters and flow fields
pex = h5read(gridName, '/pex');
k = wave_n*pex; % Wavenumber in x-direction
time = double(h5read(fname, '/time'));
a=ak/k;

fprintf('Reading %s...\n', fname);
u = h5read(fname, '/u');
v = h5read(fname, '/v');
w = h5read(fname, '/w');
eta = h5read(fname,'/eta');

% Extract the surface layer (assuming z-index 1 is the wavy boundary)
u_surf = double(squeeze(u(:,:,1)));
v_surf = double(squeeze(v(:,:,1)));
w_surf = double(squeeze(w(:,:,1)));

% Generate 2D X-coordinate grid matching the horizontal slice
Lx = 2*pi / pex;
x_1d = (0:Nx-1) * (Lx/Nx);
[X, ~] = ndgrid(x_1d, 1:Ny);

% Calculate theoretical Airy wave surface velocities from Equation 2.23
u_theory = ak * c * cos(k * (X - c*time));
v_theory = zeros(size(X));
w_theory = ak * c * sin(k * (X - c*time));
eta_theory= a* cos(k * (X - c*time));
% Compute maximum absolute errors across the entire surface
err_u = max(abs(u_surf(:) - u_theory(:)));
err_v = max(abs(v_surf(:) - v_theory(:)));
err_w = max(abs(w_surf(:) - w_theory(:)));

fprintf('\nBoundary Condition Error (Max Absolute Difference):\n');
fprintf('u-velocity error: %e\n', err_u);
fprintf('v-velocity error: %e\n', err_v);
fprintf('w-velocity error: %e\n', err_w);

% Optional: Plot a 1D slice along x (at y-index 1) to visually confirm phase alignment
figure;

subplot(1,2,1)
hold on
plot(x_1d,   u_surf(:,1), 'k-', 'LineWidth', 1.5); 
plot(x_1d, u_theory(:,1), 'r--', 'LineWidth', 1.5);
hold off
title('Streamwise Surface Velocity: Simulation vs. Eq 2.23');
xlabel('x'); ylabel('u_s');
legend('HDF5 Data', 'Theory');
grid on;

subplot(1,2,2)
hold on
plot(x_1d,   w_surf(:,1), 'k-', 'LineWidth', 1.5); 
plot(x_1d, w_theory(:,1), 'r--', 'LineWidth', 1.5);
hold off
title('Streamwise Surface Velocity: Simulation vs. Eq 2.23');
xlabel('x'); ylabel('u_s');
legend('HDF5 Data', 'Theory');
grid on;

figure
hold on
plot(x_1d,eta(:,1),'k-')
plot(x_1d,eta_theory(:,1),'--r')
hold off