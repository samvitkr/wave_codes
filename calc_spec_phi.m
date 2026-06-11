clear
close all
baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
Nx=256;
Ny=192;
Nz=128;
c=2;
fn    = 'grid.h5';
fname = fullfile(baseDir,fn);
pex  = h5read(fname, '/pex');
pey  = h5read(fname, '/pey');
kx=pex*[0:Nx/2-1,-Nx/2:-1]';
ky=pey*[0:Ny/2-1,-Ny/2:-1]';
Lx=2*pi/pex;
Ly=2*pi/pey;
%L0=Lx/wave_n;
x=[0:Nx-1]*Lx/Nx;

load(fullfile(baseDir,"slines.mat"));
eta= a0*cos(k0 * x_1D);
% sl = cos(k0 * x_1D) * asl + msl;
% sldydx = -sin(k0 * x_1D) * (asl*k0);
% Jarclength=sqrt(1+sldydx.^2);
H_bar=1;
[xi_clean, idx_x]   = unique(xi_vec);
[psi_clean, idx_y]  = unique(psi_vec);
[zeta_clean, idx_z] = unique(zeta_vec);
[XI_3D, PSI_3D, ZETA_3D] = ndgrid(xi_clean, psi_clean, zeta_clean);

Fu_3D   = griddedInterpolant(XI_3D, PSI_3D, ZETA_3D, zeros(size(XI_3D)), 'linear', 'none');

XQ = double(permute(repmat(xq, [1, 1, Ny]), [1, 3, 2]));
YQ = double(repmat(psi_clean(:)', [Nx, 1, length(msl)]));
numl=length(msl);
zq = double(slq);

ZQ = double(permute(repmat(zq, [1, 1, Ny]), [1, 3, 2]));
ETA_Q  = double(a0 * cos(k0 * XQ));
ZETA_Q = double((ZQ - ETA_Q) ./ (H_bar - ETA_Q));
%load("potvel.mat");
load(fullfile(baseDir, 'phi_interp_2d.mat'));
up=sqrt(uphi.^2+wphi.^2);
uph=uphi./up;
wph=wphi./up;
W=1./sqrt(up);

tic

tstart=3025000000;
step  =   5000000;
tend  =3820000000;

me = [];
counter=0;

phiad=zeros(Nx,Ny,numl);
phist=zeros(Nx,Ny,numl) ;


%    fn=sprintf('Sol%014d.h5',tstart);
%    fname = fullfile(baseDir,fn);
% pex = h5read(fname,'/pex');
%    pey = h5read(fname,'/pey');
%    kx=pex*[0:Nx/2-1,0,-Nx/2+1:-1]';



for tstep=tstart:step:tend
    counter=counter+1;
    time = tstep/(1e+8)
    ct=c*time;
    
    fn=sprintf('Sol%014d.h5',tstep);
    fnmat=sprintf('gradflux%014d.mat',tstep);
    fname = fullfile(baseDir,fn);
    fnamemat=fullfile(baseDir,fnmat);
    mt=matfile(fnamemat);

    u = h5read(fname, '/u')-c;
    v = h5read(fname, '/v');
    w = mt.wc;
    ox = mt.dwdy-mt.dvdz;
    oy = mt.dudz-mt.dwdx;
    oz = mt.dvdx-mt.dudy;
    clear mt
	
    kd = exp((1i*ct).*kx);
    kdis = reshape(kd,[Nx,1,1]);

    u =  ifft( (fft(u,[],1).*kdis),[],1,'symmetric');
    v =  ifft( (fft(v,[],1).*kdis),[],1,'symmetric');
    w =  ifft( (fft(w,[],1).*kdis),[],1,'symmetric');
    ox =  ifft( (fft(ox,[],1).*kdis),[],1,'symmetric');
    oy =  ifft( (fft(oy,[],1).*kdis),[],1,'symmetric');
    oz =  ifft( (fft(oz,[],1).*kdis),[],1,'symmetric');

    %ua = W.*(u.*uph + w.*wph);
    %ub = W.*v;
    %uc = W.*(-u.*wph + w.*uph);
    %oa = W.*(ox.*uph + oz.*wph);
    %ob = W.*oy;
    %oc = W.*(-ox.*wph + oz.*uph);

    %Fu_3D.Values=W.*(u.*uph + w.*wph);
    %ua =Fu_3D(XQ, YQ, ZETA_Q);
    Fu_3D.Values=W.*v;
    ub = Fu_3D(XQ, YQ, ZETA_Q);
    Fu_3D.Values= W.*(-u.*wph + w.*uph);
    uc = Fu_3D(XQ, YQ, ZETA_Q);
    %Fu_3D.Values= W.*(ox.*uph + oz.*wph);
    %oa =Fu_3D(XQ, YQ, ZETA_Q);
    Fu_3D.Values= W.*oy;
    ob = Fu_3D(XQ, YQ, ZETA_Q);
    Fu_3D.Values= W.*(-ox.*wph + oz.*uph);
    oc = Fu_3D(XQ, YQ, ZETA_Q);
    clear ox oy oz u v w

    fub=fft2(ub);
    fuc=fft2(uc);
    fob=fft2(ob);
    foc=fft2(oc);
    clear ub uc ob oc

    phiad = phiad - fuc.*(conj(fob))./(Nx*Ny);
    phist = phist + fub.*(conj(foc))./(Nx*Ny);
    clear fub fuc fob foc
end
ky=pey*[0:Ny/2-1,0,-Ny/2+1:-1]';
dkx=pex;
dky=pey;
Lx=2*pi/pex;
Ly=2*pi/pey;
dX=Lx/Nx;
dY=Ly/Ny;
% Calculate the correct continuous rescaling factor
scale_factor = (dX * dY) / (dkx * dky);


phiad=phiad./counter;
phist=phist./counter;
%%
zz   = h5read(fname, '/zz');

% Rescale the discrete co-spectra into a continuous spectral density
phiad = phiad.*scale_factor;
phist = phist.*scale_factor;
%%

%% PERFECT SPECTRAL FOLDING (For Statistical Convergence)
% Enforces symmetries to smooth the co-spectra by properly aligning
% and summing the positive and negative wavenumbers into a 1-sided spectrum.

% Pre-allocate the folded arrays (1st quadrant: kx >= 0, ky >= 0)
phiadm = zeros(Nx/2, Ny/2, numl);
phistm = zeros(Nx/2, Ny/2, numl);

% Define the exact index ranges for matching positive and negative wavenumbers
pos_x = 2:Nx/2;  neg_x = Nx:-1:(Nx/2+2);
pos_y = 2:Ny/2;  neg_y = Ny:-1:(Ny/2+2);

% --- FOLDING PHIAD ---
% 1. The (0,0) mode: No folding
phiadm(1, 1, :) = phiad(1, 1, :);

% 2. The kx axis (ky = 0): Fold +kx and -kx
phiadm(pos_x, 1, :) = phiad(pos_x, 1, :) + phiad(neg_x, 1, :);

% 3. The ky axis (kx = 0): Fold +ky and -ky
phiadm(1, pos_y, :) = phiad(1, pos_y, :) + phiad(1, neg_y, :);

% 4. The 2D Quadrants (kx > 0, ky > 0): Fold all 4 quadrants together
phiadm(pos_x, pos_y, :) = phiad(pos_x, pos_y, :) ...         % (+kx, +ky)
                        + phiad(neg_x, pos_y, :) ...         % (-kx, +ky)
                        + phiad(pos_x, neg_y, :) ...         % (+kx, -ky)
                        + phiad(neg_x, neg_y, :);            % (-kx, -ky)

% --- FOLDING PHIST ---
% 1. The (0,0) mode: No folding
phistm(1, 1, :) = phist(1, 1, :);

% 2. The kx axis (ky = 0): Fold +kx and -kx
phistm(pos_x, 1, :) = phist(pos_x, 1, :) + phist(neg_x, 1, :);

% 3. The ky axis (kx = 0): Fold +ky and -ky
phistm(1, pos_y, :) = phist(1, pos_y, :) + phist(1, neg_y, :);

% 4. The 2D Quadrants (kx > 0, ky > 0): Fold all 4 quadrants together
phistm(pos_x, pos_y, :) = phist(pos_x, pos_y, :) ...
                        + phist(neg_x, pos_y, :) ...
                        + phist(pos_x, neg_y, :) ...
                        + phist(neg_x, neg_y, :);

phinlm = phiadm + phistm;
% Strip residual imaginary floating-point noise from the FFTs
phinlm = real(phinlm);
phiadm = real(phiadm);
phistm = real(phistm);

% --- FIX 3: Truncate kx and ky to match the 1st quadrant dimensions ---
kphi = kx(1:Nx/2);
ky_folded = ky(1:Ny/2);

fphi=fullfile(baseDir,'spectra_phi.mat');
save(fphi,'phinlm','phiadm','phistm','kphi','ky_folded','zz');
disp('Co-spectra successfully folded and saved!');
