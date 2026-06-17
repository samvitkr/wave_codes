clear; close all; clc;

% =========================================================================
% 1. SETUP & GEOMETRY
% =========================================================================
ak = 0.2;
wave_n = 4; % Number of waves
pex = 0.5;
Lx = 2*pi/pex;
k0 = wave_n*pex;
a = ak/k0;
H = 1;

% Grid Resolution
Nx = 128;      % Uniform nodes in X
M  = 64;       % Chebyshev nodes in Y

xMin = 0; xMax = Lx;

% Computational Grid
xi = Lx * (0:Nx)' / Nx; % Size [Nx+1 x 1] (Uniform grid)
m_vec = (0:M)';
eta = cos(pi * m_vec / M); % Size [M+1 x 1] (Chebyshev grid)

[XI, ETA] = meshgrid(xi, eta); % Size [(M+1) x (Nx+1)]

% =========================================================================
% 2. GEOMETRIC MAPPING & METRICS
% =========================================================================
hb = a * cos(k0 * XI);
dhb_dxi = -a * k0 * sin(k0 * XI);
Y_phys = ETA .* (H - hb)/2 + (H + hb)/2;

% Mapping Metrics
eta_y = 2 ./ (H - hb);
eta_x = ((ETA - 1) ./ (H - hb)) .* dhb_dxi;

% Second-order mapping metrics
d2hb_dxi2 = -a * k0^2 * cos(k0 * XI);
detax_deta = dhb_dxi ./ (H - hb);
detax_dxi = ((ETA - 1) ./ (H - hb)) .* d2hb_dxi2 + ...
            ((ETA - 1) ./ (H - hb).^2) .* (dhb_dxi.^2);

beta_coef = detax_dxi + eta_x .* detax_deta;

% =========================================================================
% 3. EXPLICIT DCT COLLOCATION MATRICES (X) AND CHEBYSHEV (Y)
% =========================================================================
% 3A. Generate Periodic Fourier Matrices for a 2Nx Mirrored Domain
D1_2N = zeros(2*Nx, 2*Nx);
D2_2N = zeros(2*Nx, 2*Nx);
for j = 1:2*Nx
    for k = 1:2*Nx
        if j == k
            D1_2N(j,k) = 0;
            D2_2N(j,k) = - (pi^2)/(3*(2*Lx)^2) * ((2*Nx)^2 + 2);
        else
            diff_x = (j - k) * (2*Lx) / (2*Nx);
            D1_2N(j,k) = (pi/(2*Lx)) * (-1)^(j-k) * cot(pi*diff_x / (2*Lx));
            
            % CORRECTED PREFACTOR
            D2_2N(j,k) = -0.5 * (pi/Lx)^2 * (-1)^(j-k) / sin(pi*diff_x / (2*Lx))^2;
        end
    end
end

% 3B. Construct the Even-Extension Matrix (E)
E = zeros(2*Nx, Nx+1);
for i = 1:2*Nx
    if i <= Nx+1
        E(i, i) = 1;
    else
        E(i, 2*Nx - i + 2) = 1; % Mirror the domain for Cosine symmetry
    end
end

% 3C. DCT Matrices 
Dx_cos = D1_2N(1:Nx+1, :) * E;
Dxx_cos = D2_2N(1:Nx+1, :) * E;

% 3D. Standard Chebyshev Matrices for Y
Dy = zeros(M+1, M+1); cy = [2; ones(M-1,1); 2];
for i = 0:M
    for j = 0:M
        if i == j
            if i == 0, Dy(i+1,j+1) = (2*M^2+1)/6;
            elseif i == M, Dy(i+1,j+1) = -(2*M^2+1)/6;
            else, Dy(i+1,j+1) = -eta(i+1)/(2*(1-eta(i+1)^2)); end
        else, Dy(i+1,j+1) = (cy(i+1)/cy(j+1))*((-1)^(i+j))/(eta(i+1)-eta(j+1)); end
    end
end
Dy2 = Dy * Dy;

% =========================================================================
% 4. ASSEMBLE DENSE SYSTEM FOR STREAM FUNCTION
% =========================================================================
DOF = (M+1) * (Nx+1);
A = zeros(DOF, DOF); b = zeros(DOF, 1);
idx = @(m_val, i_val) (i_val) * (M+1) + (m_val + 1);

fprintf('Assembling DCT Collocation Matrix for Stream Function...\n');

for i = 0:Nx
    for m_idx = 0:M
        row = idx(m_idx, i);
        
        if m_idx == 0
            A(row, row) = 1; b(row) = H; % Top Boundary
        elseif m_idx == M
            A(row, row) = 1; b(row) = 0; % Bottom Boundary
        else
            % Interior & Side Boundaries (Neumann natively enforced by Dx_cos)
            c_cross = 2 * eta_x(m_idx+1, i+1);
            c_deta2 = eta_x(m_idx+1, i+1)^2 + eta_y(m_idx+1, i+1)^2;
            c_deta1 = beta_coef(m_idx+1, i+1);
            
            for k = 0:Nx
                A(row, idx(m_idx, k)) = A(row, idx(m_idx, k)) + Dxx_cos(i+1, k+1);
            end
            for k = 0:M
                A(row, idx(k, i)) = A(row, idx(k, i)) + c_deta2 * Dy2(m_idx+1, k+1) + c_deta1 * Dy(m_idx+1, k+1);
            end
            for kx = 0:Nx
                for ky = 0:M
                    A(row, idx(ky, kx)) = A(row, idx(ky, kx)) + c_cross * Dx_cos(i+1, kx+1) * Dy(m_idx+1, ky+1);
                end
            end
            b(row) = 0;
        end
    end
end

% =========================================================================
% 5. SOLVE & EXTRACT VELOCITIES
% =========================================================================
fprintf('Solving Stream Function System...\n');
psi_vec = A \ b;
psi = reshape(psi_vec, [M+1, Nx+1]);

dpsi_dxi = zeros(M+1, Nx+1);
for m_idx = 0:M
    dpsi_dxi(m_idx+1, :) = (Dx_cos * psi(m_idx+1, :).').';
end
dpsi_deta = Dy * psi;

u = eta_y .* dpsi_deta;
w = -(dpsi_dxi + eta_x .* dpsi_deta);

% =========================================================================
% 6. RECONSTRUCT POTENTIAL (PHI) VIA SINE TRANSFORM INTEGRATION
% =========================================================================
fprintf('Integrating U to construct exact Dirichlet Phi field...\n');

u_minus_1 = u - 1;
varphi = zeros(M+1, Nx+1);
k_vec = [0:Nx-1, 0, -Nx+1:-1]' * (pi / Lx); 

for m_idx = 0:M
    u_row = u_minus_1(m_idx+1, :);
    u_ext = [u_row, u_row(end-1:-1:2)]; % Create even extension
    
    u_hat = fft(u_ext);
    
    % Safe integration: prevent division by zero at DC and Nyquist modes
    k_safe = k_vec.';
    k_safe(k_safe == 0) = 1; 
    
    varphi_hat = u_hat ./ (1i * k_safe);
    varphi_hat(k_vec.' == 0) = 0; % Force integration constants strictly to zero
    
    varphi_ext = real(ifft(varphi_hat));
    varphi(m_idx+1, :) = varphi_ext(1:Nx+1); % Extracts exact Sine series
end

phi = varphi + XI;

% =========================================================================
% 7. ALIAS-FREE DIVERGENCE CHECK
% =========================================================================
du_dxi = zeros(M+1, Nx+1);
for m_idx = 0:M
    du_dxi(m_idx+1, :) = (Dx_cos * u(m_idx+1, :).').'; 
end

du_deta = Dy * u;
dw_deta = Dy * w;

du_dx = du_dxi + eta_x .* du_deta;
dw_dy = eta_y .* dw_deta;
divergence_field = du_dx + dw_dy;

max_div = max(abs(divergence_field(:)));

fprintf('\n==================================================\n');
fprintf('  SYSTEM VERIFICATION (DCT METHOD)\n');
fprintf('==================================================\n');
fprintf('  Max Divergence Error:       %e\n', max_div);
fprintf('  Max Inlet Phi Fluctuation:  %e\n', max(abs(phi(:,1) - 0)));
fprintf('  Max Outlet Phi Fluctuation: %e\n', max(abs(phi(:,end) - Lx)));
fprintf('==================================================\n\n');

% =========================================================================
% 8. VISUALIZATION
% =========================================================================
%figure('Color', 'w', 'Position', [100, 100, 1200, 450]);

%% subplot(2,1,1);
figure
contourf(XI, Y_phys, phi, 40, 'EdgeColor', 'none'); hold on;
plot(XI(end, :), Y_phys(end, :), 'k', 'LineWidth', 2); 
plot(XI(1, :), Y_phys(1, :), 'k', 'LineWidth', 2);   
colormap('jet'); colorbar; axis equal;
xlim([xMin xMax]); ylim([-a-0.05 H+0.05]);
title('Velocity Potential (\phi): Exact Dirichlet Boundaries');
xlabel('x'); ylabel('z');

%subplot(2,1,2);
figure
contourf(XI, Y_phys, psi, 40, 'EdgeColor', 'none'); hold on;
plot(XI(end, :), Y_phys(end, :), 'k', 'LineWidth', 2); 
plot(XI(1, :), Y_phys(1, :), 'k', 'LineWidth', 2);   
colormap('jet'); colorbar; axis equal;
xlim([xMin xMax]); ylim([-a-0.05 H+0.05]);
title('Stream Function (\psi): Zero Corner Singularity');
xlabel('x'); ylabel('z');


