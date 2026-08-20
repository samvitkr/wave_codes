clear;
close all;

% ==========================================
% 1. PARAMETERS & INITIALIZATION
% ==========================================
baseDir = '/scratch.global/kuma0458/c8ak1_re180/run';
c       = 8;
ret     = 180;
nu      = 1 / ret;

Nx      = 256;
Ny      = 192;
Nz      = 128;

tstart  = 4021250000; 
step    = 1250000;
tend    = 4270000000;

% ==========================================
% 2. LOAD STATIC GRID & METRICS
% ==========================================
fname_gridh5 = fullfile(baseDir, 'grid.h5');
fprintf('Loading grid definitions...\n');
zz  = h5read(fname_gridh5, '/zz');
zw  = h5read(fname_gridh5, '/zw');
pex = h5read(fname_gridh5, '/pex');
pey = h5read(fname_gridh5, '/pey');

Lx = 2 * pi / pex;
Ly = 2 * pi / pey;

% Wavenumbers (reshaped for implicit expansion)
kx    = pex * [0:Nx/2-1, -Nx/2:-1]';
ky    = pey * [0:Ny/2-1, -Ny/2:-1]';
kx_3D = reshape(kx, [Nx, 1, 1]);
ky_3D = reshape(ky, [1, Ny, 1]);

% Load static grid file (Since wave is frozen, t=0 metrics are exact)
load(fullfile(baseDir, 'grid.mat'), 'X', 'Y', 'Z', 'Zw', 'dZetadz', 'dZetadx');

dx = X(4,4,4) - X(3,3,3);
dy = Y(4,4,4) - Y(3,3,3);

% Ensure proper 3D arrays for curvilinear metric terms
dZetadx_3D = reshape(dZetadx, [Nx, 1, Nz]);
dZetadz_3D = reshape(dZetadz, [Nx, 1, 1]); 
Jacobian   = 1 ./ dZetadz_3D;
Jac_1D     = squeeze(Jacobian(:, 1, 1)); 
vol        = trapz(X(:, 1, 1), Jac_1D);

% ==========================================
% 3. LOAD AUXILIARY DATA
% ==========================================
load(fullfile(baseDir, 'flowrate.mat'), 't', 'flowrate', 'Jdot');
load(fullfile(baseDir, 'phi_interp_2d.mat'), 'uphi', 'wphi', 'J');

uphi = uphi ./ J;
wphi = wphi ./ J;

Ts      = Jdot .* 0;
Tnls    = Jdot .* 0;
Tviscs  = Jdot .* 0;
Tconvs  = Jdot .* 0;
Tstrs   = Jdot .* 0;
phidots = Jdot .* 0;
check   = Jdot .* 0;

% Z-Derivative Stencil Setup
dz = diff(zz);
nz = Nz;
idx = 1:(nz-2);
dz_i   = dz(idx);
dz_ip1 = dz(idx+1);
sumdz  = dz_i + dz_ip1;

a_st = reshape(1./dz_ip1 - 1./sumdz, 1, 1, []);
b_st = reshape(1./dz_i   - 1./dz_ip1, 1, 1, []);
c_st = reshape(1./sumdz  - 1./dz_i, 1, 1, []);

% ==========================================
% 4. MAIN TIME LOOP
% ==========================================
for tstep = tstart:step:tend
    fname = fullfile(baseDir, sprintf('Sol%014d.h5', tstep));
    fprintf('Processing %s\n', fname);
    
    u    = h5read(fname, '/u');
    v    = h5read(fname, '/v');
    w    = h5read(fname, '/w');
    time = h5read(fname, '/time');
    ct   = c * time;
    
    % --- A. SHIFT VELOCITIES TO STATIONARY WAVE FRAME FIRST ---
    kdis = reshape(exp((1i * ct) .* kx), [Nx, 1, 1]);
    
    u = ifft(fft(u, [], 1) .* kdis, [], 1, 'symmetric') - c;
    v = ifft(fft(v, [], 1) .* kdis, [], 1, 'symmetric');
    w = ifft(fft(w, [], 1) .* kdis, [], 1, 'symmetric');
    
    % Interpolate w to cell centers
    wc = permute(w, [3 1 2]);
    wc = interp1(zw(1:end-1), wc(1:end-1, :, :), zz);
    wc = single(permute(wc, [2 3 1]));
    
    % --- B. VERTICAL GRADIENTS (d/dzeta) ---
    dudzeta = zeros(size(u), 'single');
    dvdzeta = zeros(size(v), 'single');
    dwdzeta = zeros(size(wc), 'single');
    
    dudzeta(:, :, 1) = (u(:, :, 2)  - u(:, :, 1))  ./ dz(1);
    dvdzeta(:, :, 1) = (v(:, :, 2)  - v(:, :, 1))  ./ dz(1);
    dwdzeta(:, :, 1) = (wc(:, :, 2) - wc(:, :, 1)) ./ dz(1);
    
    dudzeta(:, :, end) = (u(:, :, end) - u(:, :, end-1)) ./ dz(end);
    dvdzeta(:, :, end) = (v(:, :, end) - v(:, :, end-1)) ./ dz(end);
    dwdzeta(:, :, end) = (wc(:, :, end)- wc(:, :, end-1))./ dz(end);
    
    dudzeta(:, :, 2:nz-1) = a_st.*u(:, :, 3:nz)  + b_st.*u(:, :, 2:nz-1)  + c_st.*u(:, :, 1:nz-2);
    dvdzeta(:, :, 2:nz-1) = a_st.*v(:, :, 3:nz)  + b_st.*v(:, :, 2:nz-1)  + c_st.*v(:, :, 1:nz-2);
    dwdzeta(:, :, 2:nz-1) = a_st.*wc(:, :, 3:nz) + b_st.*wc(:, :, 2:nz-1) + c_st.*wc(:, :, 1:nz-2);
    
    dudz = single(dZetadz_3D .* dudzeta);
    dvdz = single(dZetadz_3D .* dvdzeta);
    dwdz = single(dZetadz_3D .* dwdzeta);
    
    % --- C. HORIZONTAL FOURIER DERIVATIVES ---
    dudy = single(ifft(fft(u, [], 2)  .* (1i .* ky_3D), [], 2, 'symmetric'));
    dvdy = single(ifft(fft(v, [], 2)  .* (1i .* ky_3D), [], 2, 'symmetric'));
    dwdy = single(ifft(fft(wc, [], 2) .* (1i .* ky_3D), [], 2, 'symmetric'));
    
    d2udy2 = single(ifft(fft(u, [], 2)  .* (-(ky_3D).^2), [], 2, 'symmetric'));
    d2vdy2 = single(ifft(fft(v, [], 2)  .* (-(ky_3D).^2), [], 2, 'symmetric'));
    d2wdy2 = single(ifft(fft(wc, [], 2) .* (-(ky_3D).^2), [], 2, 'symmetric'));
    
    dudxi = ifft(fft(u, [], 1)  .* (1i .* kx_3D), [], 1, 'symmetric');
    dvdxi = ifft(fft(v, [], 1)  .* (1i .* kx_3D), [], 1, 'symmetric');
    dwdxi = ifft(fft(wc, [], 1) .* (1i .* kx_3D), [], 1, 'symmetric');
    
    dudx = single(dudxi + dZetadx_3D .* dudzeta);
    dvdx = single(dvdxi + dZetadx_3D .* dvdzeta);
    dwdx = single(dwdxi + dZetadx_3D .* dwdzeta);
    
    % --- D. EXACT VISCOUS LAPLACIANS ---
    Jinv = 1 ./ dZetadz_3D;
    
    F1_u = Jinv .* dudx;
    dF1_u_dxi = ifft(fft(F1_u, [], 1) .* (1i .* kx_3D), [], 1, 'symmetric');
    F3_u = Jinv .* (dZetadx_3D .* dudx + dZetadz_3D.^2 .* dudzeta);
    dF3_u_dzeta = zeros(size(u), 'single');
    dF3_u_dzeta(:, :, 1)     = (F3_u(:, :, 2) - F3_u(:, :, 1)) ./ dz(1);
    dF3_u_dzeta(:, :, 2:nz-1) = a_st.*F3_u(:, :, 3:nz) + b_st.*F3_u(:, :, 2:nz-1) + c_st.*F3_u(:, :, 1:nz-2);
    dF3_u_dzeta(:, :, end)   = (F3_u(:, :, end) - F3_u(:, :, end-1)) ./ dz(end);
    viscu = single(nu .* ((dF1_u_dxi + dF3_u_dzeta) ./ Jinv + d2udy2));
    
    F1_w = Jinv .* dwdx;
    dF1_w_dxi = ifft(fft(F1_w, [], 1) .* (1i .* kx_3D), [], 1, 'symmetric');
    F3_w = Jinv .* (dZetadx_3D .* dwdx + dZetadz_3D.^2 .* dwdzeta);
    dF3_w_dzeta = zeros(size(wc), 'single');
    dF3_w_dzeta(:, :, 1)     = (F3_w(:, :, 2) - F3_w(:, :, 1)) ./ dz(1);
    dF3_w_dzeta(:, :, 2:nz-1) = a_st.*F3_w(:, :, 3:nz) + b_st.*F3_w(:, :, 2:nz-1) + c_st.*F3_w(:, :, 1:nz-2);
    dF3_w_dzeta(:, :, end)   = (F3_w(:, :, end) - F3_w(:, :, end-1)) ./ dz(end);
    viscw = single(nu .* ((dF1_w_dxi + dF3_w_dzeta) ./ Jinv + d2wdy2));
    
    % --- E. VORTICITY & JA COMPUTATIONS ---
    ox = dwdy - dvdz;
    oy = dudz - dwdx;
    oz = dvdx - dudy;
    
    voz = single(v .* oz);
    woy = single(wc .* oy);
    uoy = single(u .* oy);  % u is already u_lab - c
    vox = single(v .* ox);
    
    JAnl   = uphi.*(voz - woy) + wphi.*(uoy - vox);
    JAconv = uphi.*(   - woy)  + wphi.*(uoy);
    JAstr  = uphi.*(voz)       + wphi.*(   - vox);
    JAvisc = uphi.*viscu       + wphi.*viscw;
    
    JAnl(isnan(JAnl))     = 0;
    JAvisc(isnan(JAvisc)) = 0;
    JAconv(isnan(JAconv)) = 0;
    JAstr(isnan(JAstr))   = 0;
    JAtot = JAnl + JAvisc;
    
    % Save JA Fields
    fn_jafields = fullfile(baseDir, sprintf('jafieldstat%014d.mat', tstep));
    save(fn_jafields, 'JAnl', 'JAvisc', 'JAconv', 'JAstr');
    
    % --- F. VOLUME INTEGRATIONS ---
    JAin     = trapz(zz, squeeze(dy*trapz(JAtot, 2))./Ly, 2);
    JAnlin   = trapz(zz, squeeze(dy*trapz(JAnl, 2))./Ly, 2);
    JAviscin = trapz(zz, squeeze(dy*trapz(JAvisc, 2))./Ly, 2);
    JAconvin = trapz(zz, squeeze(dy*trapz(JAconv, 2))./Ly, 2);
    JAstrin  = trapz(zz, squeeze(dy*trapz(JAstr, 2))./Ly, 2);
    
    T     = -trapz(X(:, 1, 1), Jac_1D .* JAin) / J;
    Tnl   = -trapz(X(:, 1, 1), Jac_1D .* JAnlin) / J;
    Tvisc = -trapz(X(:, 1, 1), Jac_1D .* JAviscin) / J;
    Tconv = -trapz(X(:, 1, 1), Jac_1D .* JAconvin) / J;
    Tstr  = -trapz(X(:, 1, 1), Jac_1D .* JAstrin) / J;
    
    it = find(abs(t - time) < 1e-6, 1);
    if ~isempty(it)
        phidot = Jdot(it) * Lx / J;
        
        Ts(it)      = T;
        Tnls(it)    = Tnl;
        Tviscs(it)  = Tvisc;
        Tconvs(it)  = Tconv;
        Tstrs(it)   = Tstr;
        phidots(it) = phidot;
        
        check(it)   = 100 * (1 - (T + phidot) / vol);
    end
end

% ==========================================
% 5. FINALIZE & SAVE
% ==========================================
mja = fullfile(baseDir, 'JAseries_static.mat');
save(mja, 'Ts', 'Tnls', 'Tviscs', 'Tconvs', 'Tstrs', 'phidots', 'check');
fprintf('Unified JA relation script completed and saved.\n');




