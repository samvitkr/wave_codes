close all
clear
baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
Nx=256;
Ny=192;
Nz=128;
wave_n=12;
ak=0.2;

Nxl = round(Nx./wave_n)+1;
load(fullfile( baseDir,'grid.mat'))
load(fullfile(baseDir,'slines.mat'));
a0=ak/k0;
x_1D = X(:,1,1);
%load(fullfile(baseDir,'mean_phase_fields.mat'))
eta= a0*cos(k0 * x_1D);
% sl = cos(k0 * x_1D) * asl + msl;
% sldydx = -sin(k0 * x_1D) * (asl*k0);
% Jarclength=sqrt(1+sldydx.^2);
H_bar=1;
%XI_2D = squeeze(X(:,1,:));
%ZETA_2D = Zsl - eta;
% [xi_clean, idx_x]   = unique(xi_vec);
% [zeta_clean, idx_z] = unique(zeta_vec);
xi_clean=squeeze(double(X(:,5,5)));
Zs = double(squeeze(Z(:,1,:)));
XQ=squeeze(double(X(:,5,:)));

ETA_Q  = double(a0 * cos(k0 * (XQ)));
ZETA_2D = double((Zs - ETA_Q) ./ (H_bar - ETA_Q));
zeta_clean= double(squeeze(ZETA_2D(1,:)));



 [XI_2D, ZETA_2D] = ndgrid(xi_clean, zeta_clean);

Fu_2D=griddedInterpolant(XI_2D,ZETA_2D,zeros(size(XI_2D)),'makima','nearest');
ZQ = Zsl;
[~,numl]=size(Zsl)
% numl=length(msl);
% ZQ = double(sl);
ZETA_Q = double((ZQ - ETA_Q) ./ (H_bar - ETA_Q));

load(fullfile(baseDir,'phi_interp_2d.mat'),'uphi','wphi','Phi_2D')
upin = squeeze(uphi);

Fu_2D.Values=upin;
up_sl = Fu_2D(XQ,ZETA_Q);

Fu_2D.Values=Phi_2D;
phi_sl = Fu_2D(XQ,ZETA_Q);

clear uphi wphi Phi_2D

load(fullfile(baseDir,'mean_fields.mat'),'viscav','strav','convav');

Fu_2D.Values=squeeze(mean(viscav,2));
visc_sl=Fu_2D(XQ,ZETA_Q)./up_sl;

Fu_2D.Values=squeeze(mean(strav,2));
str_sl =Fu_2D(XQ,ZETA_Q)./up_sl;

Fu_2D.Values=squeeze(mean(convav,2));
conv_sl=Fu_2D(XQ,ZETA_Q)./up_sl;


 convmf=zeros(numl,1);
 viscmf = zeros(numl,1);
 strmf= zeros(numl,1);
 
 for k =1:numl
 
 	%convm(k) = trapz(phi_sl(:,k),conv_sl(:,k));
 	%strm(k) = trapz(phi_sl(:,k),str_sl(:,k));
 	%viscm(k) = trapz(phi_sl(:,k),visc_sl(:,k));
 
 	convmf(k) = trapz(XQ(:,k),conv_sl(:,k));
         strmf(k) = trapz(XQ(:,k),str_sl(:,k));
         viscmf(k) = trapz(XQ(:,k),visc_sl(:,k));
 end
% 
% %%
% 
% close all
% fname=fullfile(baseDir,'grid.h5');
% pex  = h5read(fname, '/pex');
% L = 2*pi/(pex);
% zqm = ZETA_Q(1,:)';
% figure
% hold on
% plot(-convm./L,zqm,'--b','LineWidth',1)
% plot(-strm./L,zqm,':b','LineWidth',1)
% plot((-convm-strm)./L,zqm,'-b','LineWidth',1)
% plot((-viscm)./L,zqm,'-r','LineWidth',1)
% plot((-strm-convm-viscm)./L,zqm,'-k','LineWidth',1)
% 
% xline(0)




%%

% 1. Truncate query trajectories to exactly one wave phase
XQ_p     = XQ(1:Nxl,:);
ZETA_Q_p = ZETA_Q(1:Nxl,:);

% 2. Create the dedicated 2D Phase-Averaged Interpolator with proper boundary clamping
[XI_2D_phase, ZETA_2D_phase] = ndgrid(xi_clean(1:Nxl), zeta_clean);
Fu_2D_phase = griddedInterpolant(XI_2D_phase, ZETA_2D_phase, zeros(size(XI_2D_phase)), 'makima', 'nearest');

% 3. Process the horizontal potential velocity denominator (u_phi^*)
load(fullfile(baseDir,'phi_interp_2d.mat'),'uphi')
uphi_phase = squeeze(uphi(1:Nxl, :)); % Truncate grid data to phase length

Fu_2D_phase.Values = uphi_phase;
uphi_sl = Fu_2D_phase(XQ_p, ZETA_Q_p);

clear uphi wphi Phi_2D

% 4. Load Phase-Averaged Fields and Apply Analytical Denominator Cancellation
load(fullfile(baseDir,'mean_phase_fields.mat'),'viscpav','strpav','convpav');

Fu_2D_phase.Values = viscpav;
visc_sl = Fu_2D_phase(XQ_p, ZETA_Q_p) ./ uphi_sl;

Fu_2D_phase.Values = strpav;
str_sl  = Fu_2D_phase(XQ_p, ZETA_Q_p) ./ uphi_sl;

Fu_2D_phase.Values = convpav;
conv_sl = Fu_2D_phase(XQ_p, ZETA_Q_p) ./ uphi_sl;

convm = zeros(numl,1);
viscm = zeros(numl,1);
strm  = zeros(numl,1);
zqm = ZETA_Q_p(1,:)';
% 5. Integrate strictly over the uniform 1D phase coordinate (Skipping k=1 wall boundary)
for k = 2:numl
    convm(k) = trapz(XQ_p(:,k), conv_sl(:,k));
    strm(k)  = trapz(XQ_p(:,k), str_sl(:,k));
    viscm(k) = trapz(XQ_p(:,k), visc_sl(:,k));
end

%% Plotting
close all
fname = fullfile(baseDir,'grid.h5');
pex   = h5read(fname, '/pex');
Lx    = 2*pi / pex;

% Normalization length must be one wave phase, not the full domain
lambda = Lx / wave_n;

zqm = ZETA_Q_p(1,:)';
figure
hold on
plot(-convm(2:end)./lambda, zqm(2:end), '--b', 'LineWidth', 1)
plot(-strm(2:end)./lambda,  zqm(2:end), ':b', 'LineWidth', 1)
plot((-convm(2:end)-strm(2:end))./lambda, zqm(2:end), '-b', 'LineWidth', 1)
plot((-viscm(2:end))./lambda, zqm(2:end), '-r', 'LineWidth', 1)
plot((-strm(2:end)-convm(2:end)-viscm(2:end))./lambda, zqm(2:end), '-k', 'LineWidth', 1)

%%
figure
hold on
plot(-convmf(2:end)./Lx, zqm(2:end), '--ob', 'LineWidth', 1)
plot(-strmf(2:end)./Lx,  zqm(2:end), ':ob', 'LineWidth', 1)
plot((-convmf(2:end)-strmf(2:end))./Lx, zqm(2:end), '-ob', 'LineWidth', 1)
plot((-viscmf(2:end))./Lx, zqm(2:end), '-or', 'LineWidth', 1)
plot((-strmf(2:end)-convmf(2:end)-viscmf(2:end))./Lx, zqm(2:end), '-ok', 'LineWidth', 1)
% 
% xline(0)
save( fullfile(baseDir, 'mean_sl.mat'),'convmf','strmf','viscmf','zqm','-v7.3'  )

%save( fullfile(baseDir, 'mean_sl.mat'),'convm','strm','viscm','zqm','-v7.3'  )
