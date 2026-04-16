clear
close all

% baseDir = fullfile('/projects', 'standard', 'shenl','naras062','JA-Samvit','Retau-10','RUN-c0');
 baseDir = '/users/1/kuma0458/wave/wavy_wall';
% baseDir = '/users/1/kuma0458/wave/wave_c_2';
%baseDir = '/users/1/kuma0458/wave/wavy_ret180';

tic
ret = 10;
nu = 1/ret;
tstart = 20200000;
step = 200000;
tend = 32000000;

fj = fullfile(baseDir, 'flowrate.mat');
load(fj);

Ts = Jdot .* 0;
phidots = Jdot .* 0;
Tnls = Ts;
Tviscs = Ts;
check = Jdot .* 0;

Nx = 128; Ny = 16; Nz = 128;
c = 0;

% --- MEMORY FIX 1: Load static files OUTSIDE the loop ---
fng = fullfile(baseDir, 'grid.mat');
load(fng); % Loads X, Z, dZetadz

fnp = fullfile(baseDir, 'potvel.mat');
load(fnp); % Loads uphi, wphi

% Cast to single immediately to halve memory footprint
uphis = single(uphi);
wphis = single(wphi);
clear uphi wphi % Remove double precision copies
it=0;
for tstep = tstart:step:tend
    it=it+1;
    fn = sprintf('Sol%014d.h5', tstep);
    fnmat = sprintf('gradflux%014d.mat', tstep);
    fname = fullfile(baseDir, fn);
    fnamemat = fullfile(baseDir, fnmat);
    fnja = fullfile(baseDir, sprintf('jafields%014d.mat', tstep));

    fprintf('Reading %s\n', fname);
    
    u = single(h5read(fname, '/u'));
    v = single(h5read(fname, '/v'));
    % w = h5read(fname, '/w'); <--- REMOVED: 'w' is never used in math below, saving ~1/3 memory
    
    time = h5read(fname, '/time');
    zz = h5read(fname, '/zz');
    pey = h5read(fname, '/pey');
    pex = h5read(fname, '/pex');
    zz = zz';
    Ly = 2*pi/pey;
    Lx = 2*pi/pex;

    % --- MEMORY FIX 2: Load variables into a struct to clear them aggressively ---
    GF = load(fnamemat);

    % Compute intermediate values and clear the ingredients immediately
    oy = single(GF.dudz - GF.dwdx);
    woy = single(GF.wc) .* oy;
    uoy = (u - c) .* oy;
    clear oy; 
    
    ox = single(GF.dwdy - GF.dvdz);
    vox = v .* ox;
    clear ox;
    
    oz = single(GF.dvdx - GF.dudy);
    voz = v .* oz;
    clear oz v u; % We no longer need u, v, or oz

    % Compute Nonlinear and Viscous Terms
    JAnl = uphis .* (voz - woy) + wphis .* (uoy - vox);
    clear voz woy uoy vox % Clear massive intermediates
    
    JAvisc = uphis .* single(GF.viscu) + wphis .* single(GF.viscw);
    clear GF % Free the rest of the gradflux structure entirely
    
    % Clean up NaNs
    JAnl(isnan(JAnl)) = 0;
    JAvisc(isnan(JAvisc)) = 0;
    
    save(fnja, 'JAnl', 'JAvisc');

    % --- MEMORY FIX 3: Skip creating JAx and JAz entirely ---
    % JAin is mathematically equal to JAnl + JAvisc.
    JAin = JAnl + JAvisc; 

    % Perform spatial averaging
    JAin = trapz(zz', squeeze(mean(JAin, 2)), 2);
    JAnlin = trapz(zz', squeeze(mean(JAnl, 2)), 2);
    JAviscin = trapz(zz', squeeze(mean(JAvisc, 2)), 2);
    
    clear JAnl JAvisc % Massive memory dump before integration
    
    Jacobian = 1 ./ dZetadz;
    T = -trapz(X(:,1,1), Jacobian .* JAin);
    Tnl = -trapz(X(:,1,1), Jacobian .* JAnlin);
    Tvisc = -trapz(X(:,1,1), Jacobian .* JAviscin);
phidot = (Jdot(it))*Lx/Ly;
Ts(it)=T;
Tnls(it)=Tnl;
Tviscs(it)=Tvisc;
phidots(it)=phidot;

    % Clean up remaining loop variables before the next step
    %clear JAin JAnlin JAviscin Jacobian T Tnl Tvisc
    
    % (Your commented out logic for phidot, Ts, Tnls, etc. remains untouched)
end

%save('JAseries.mat', 'Ts', 'Tnls', 'Tviscs', 'phidots', 'check')
