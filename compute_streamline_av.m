clear
close all
baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';

load("slines.mat");
%load("potvel.mat")
eta= a0*cos(k0 * x_1D);
sl = cos(k0 * x_1D) * asl + msl;
sldydx = -sin(k0 * x_1D) * (asl*k0);
Jarclength=sqrt(1+sldydx.^2);
H_bar=1;

[xi_clean, idx_x]   = unique(xi_vec);
[psi_clean, idx_y]  = unique(psi_vec);
[zeta_clean, idx_z] = unique(zeta_vec);

    [XI_3D, PSI_3D, ZETA_3D] = ndgrid(xi_clean, psi_clean, zeta_clean);

Fnl_3D   = griddedInterpolant(XI_3D, PSI_3D, ZETA_3D, zeros(size(XI_3D)), 'linear', 'none');
Fvisc_3D = griddedInterpolant(XI_3D, PSI_3D, ZETA_3D, zeros(size(XI_3D)), 'linear', 'none');
upot = griddedInterpolant(XI_3D, PSI_3D, ZETA_3D, zeros(size(XI_3D)), 'linear', 'none');

dx = xi_vec(2)-xi_vec(1);
dy = psi_vec(2)-psi_vec(1);
Nx = length(xi_vec);
Ny = length(psi_clean);
XQ = permute(repmat(xq, [1, 1, Ny]), [1, 3, 2]);
YQ = repmat(psi_clean(:)', [length(xq), 1, length(msl)]);
numl=length(msl);
nlav=zeros(numl,1);
viscav=zeros(numl,1);

%baseDir = '/users/1/kuma0458/wave/wavy_wall';
%baseDir = '/users/1/kuma0458/wave/wavy_ret180';
%baseDir = '/users/1/kuma0458/wave/wave_c_2';
%baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';

% tstart=30000000;
% step=    2000000;
% tend=  32000000;
tstart=20200000;
step=  200000;
tend=  32000000;
counter=0;
for tstep=tstart:step:tend
    % fn=               sprintf('DAT000%03d99999999',tstep)
    % fnmat= sprintf('gradflux000%03d99999999.mat',tstep)
    fngs = sprintf('jafields%014d.mat',tstep)
    fng = fullfile(baseDir,fngs);
    fnp=sprintf('potvel%014d.mat',tstep);
%fnp=sprintf('potvel.mat');
    fnp=fullfile(baseDir,fnp);
    mp=matfile(fnp);
    up=double(repmat(mp.uphi , 1,Ny,1));

    % fngrid=sprintf('grid%014d.mat',tstep)
    % fngr = fullfile(baseDir,fngrid);
    % load(fngr)
    load(fng);
    %fng=fullfile(baseDir,fng);
    t=tstep/10^8;
    
    sldis = cos(k0 * (x_1D - c*t)) * asl + msl;
    zq = double(sldis);
    ZQ = permute(repmat(zq, [1, 1, Ny]), [1, 3, 2]);
    sldydx = -sin(k0 * (x_1D-c*t)) * (asl*k0);
    Jarclen=sqrt(1+sldydx.^2);
    Jl = reshape(Jarclen,[Nx 1 length(msl)]);

    JAnl_clean   = double(JAnl(idx_x, idx_y, idx_z));
    JAvisc_clean = double(JAvisc(idx_x, idx_y, idx_z));
    
    % --- B. Update Interpolant Values (NO REBUILDING GRID) ---
    Fnl_3D.Values   = JAnl_clean;
    Fvisc_3D.Values = JAvisc_clean;
    upot.Values=up;
    %upot.Values=up;
    % --- C. Shift the Streamline Sheets ---
    
    % --- D. Map to Computational Space ---
    ETA_Q  = a0 * cos(k0 * (XQ - c*t));
    ZETA_Q = (ZQ - ETA_Q) ./ (H_bar - ETA_Q);
    
    % --- E. Instantly Extract 3D Data ---
    janlsl_3D   = Fnl_3D(XQ, YQ, ZETA_Q);
    javiscsl_3D = Fvisc_3D(XQ, YQ, ZETA_Q);
    uphi=upot(XQ,YQ,ZETA_Q);

    nlav = nlav+(dx*dy)*squeeze(sum(janlsl_3D./uphi,[1 2]));
    viscav = viscav+(dx*dy)*squeeze(sum(javiscsl_3D./uphi,[1 2]));
counter =counter+1;
    %total_J_3D  = janlsl_3D + javiscsl_3D;
    clear data JAnl_clean JAvisc_clean janlsl_3D javiscsl_3D;
end
Lx=2*pi;
Lz=0.5*pi;
nlav=nlav./(Lx*Lz*counter);
viscav=viscav./(Lx*Lz*counter);
nlav=nlav./uphi_in;
viscav=viscav./uphi_in;
