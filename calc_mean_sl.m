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
upin = squeeze(uphi.^2 + wphi.^2);

Fu_2D.Values=upin;
up_sl = Fu_2D(XQ,ZETA_Q);


clear uphi wphi Phi_2D

load(fullfile(baseDir,'mean_fields.mat'),'viscav','strav','convav','uav','oyav');

Fu_2D.Values=squeeze(mean(viscav,2));
visc_sl=Fu_2D(XQ,ZETA_Q)./up_sl;

Fu_2D.Values=squeeze(mean(strav,2));
str_sl =Fu_2D(XQ,ZETA_Q)./up_sl;

Fu_2D.Values=squeeze(mean(convav,2));
conv_sl=Fu_2D(XQ,ZETA_Q)./up_sl;

Fu_2D.Values=squeeze(mean(uav,2));
u_sl=Fu_2D(XQ,ZETA_Q)./up_sl;

Fu_2D.Values=squeeze(mean(oyav,2));
oy_sl=Fu_2D(XQ,ZETA_Q)./up_sl;

 convmf=zeros(numl,1);
 viscmf = zeros(numl,1);
 strmf= zeros(numl,1);
 umf=zeros(numl,1);
 oymf = zeros(numl,1);

 for k =1:numl
 
 	convmf(k)= trapz(phisl(:,k),conv_sl(:,k));
 	strmf(k) = trapz(phisl(:,k), str_sl(:,k));
 	viscmf(k)= trapz(phisl(:,k),visc_sl(:,k));
  	umf(k)   = trapz(phisl(:,k),   u_sl(:,k));
 	oymf(k)  = trapz(phisl(:,k),  oy_sl(:,k));

 	% convmf(k) = trapz(XQ(:,k),conv_sl(:,k));
    %      strmf(k) = trapz(XQ(:,k),str_sl(:,k));
    %      viscmf(k) = trapz(XQ(:,k),visc_sl(:,k));
 end

 zqm = ZETA_Q(1,:)';



%% Plotting
close all
fname = fullfile(baseDir,'grid.h5');
pex   = h5read(fname, '/pex');
Lx    = 2*pi / pex;
umf =   umf./Lx;
oymf = oymf./Lx;
% Normalization length must be one wave phase, not the full domain
lambda = Lx / wave_n;

%zqm = ZETA_Q_p(1,:)';
% figure
% hold on
% plot(-convm(2:end)./lambda, zqm(2:end), '--b', 'LineWidth', 1)
% plot(-strm(2:end)./lambda,  zqm(2:end), ':b', 'LineWidth', 1)
% plot((-convm(2:end)-strm(2:end))./lambda, zqm(2:end), '-b', 'LineWidth', 1)
% plot((-viscm(2:end))./lambda, zqm(2:end), '-r', 'LineWidth', 1)
% plot((-strm(2:end)-convm(2:end)-viscm(2:end))./lambda, zqm(2:end), '-k', 'LineWidth', 1)

%%
figure
hold on
plot(-convmf(2:end)./Lx, zqm(2:end), '--b', 'LineWidth', 1)
plot(-strmf(2:end)./Lx,  zqm(2:end), ':b', 'LineWidth', 1)
plot((-convmf(2:end)-strmf(2:end))./Lx, zqm(2:end), '-b', 'LineWidth', 1)
plot((-viscmf(2:end))./Lx, zqm(2:end), '-r', 'LineWidth', 1)
plot((-strmf(2:end)-convmf(2:end)-viscmf(2:end))./Lx, zqm(2:end), '-k', 'LineWidth', 1)
% 
% xline(0)
save( fullfile(baseDir, 'mean_sl.mat'),'convmf','strmf','viscmf','umf','oymf','zqm','-v7.3'  )

%save( fullfile(baseDir, 'mean_sl.mat'),'convm','strm','viscm','zqm','-v7.3'  )
%%
figure
subplot(1,2,1)
plot(zqm,umf)
subplot(1,2,2)
plot(zqm,oymf)

% load(fullfile(baseDir,'mean_phase_fields.mat'));
% x_1D = Xph(:,1);
% eta= a0*cos(k0 * x_1D);
% H_bar=1;
% xi_clean=x_1D;
% Zs = double(squeeze(Zph(:,:)));
% XQ = squeeze(double(Xph(:,:)));
% 
% ETA_Q  = double(a0 * cos(k0 * (XQ)));
% ZETA_2D = double((Zs - ETA_Q) ./ (H_bar - ETA_Q));
% zeta_clean= double(squeeze(ZETA_2D(1,:)));
% [XI_2D, ZETA_2D] = ndgrid(xi_clean, zeta_clean);
% 
% Fu_2D=griddedInterpolant(double(XI_2D),double(ZETA_2D),ZETA_2D.*0 ,'makima','nearest');
% ZQ = Zslph;
% [~,numl]=size(Zslph);
% ZETA_Q = double((ZQ - ETA_Q) ./ (H_bar - ETA_Q));
% upin = squeeze(uphiph.^2 + wphiph.^2);
% 
% Fu_2D.Values=upin;
% up_sl = Fu_2D(XQ,ZETA_Q);
% 
% Fu_2D.Values=viscph;
% visc_sl=Fu_2D(XQ,ZETA_Q)./up_sl;
% 
% Fu_2D.Values=strph;
% str_sl =Fu_2D(XQ,ZETA_Q)./up_sl;
% 
% Fu_2D.Values=convph;
% conv_sl=Fu_2D(XQ,ZETA_Q)./up_sl;
% 
% convmfp=zeros(numl,1);
%  viscmfp = zeros(numl,1);
%  strmfp= zeros(numl,1);
% 
%  for k =1:numl
% 
%  	convmfp(k) = trapz(phislph(:,k),conv_sl(:,k));
%  	strmfp(k) = trapz(phislph(:,k),str_sl(:,k));
%  	viscmfp(k) = trapz(phislph(:,k),visc_sl(:,k));
% 
%  	% convmf(k) = trapz(XQ(:,k),conv_sl(:,k));
%     %      strmf(k) = trapz(XQ(:,k),str_sl(:,k));
%     %      viscmf(k) = trapz(XQ(:,k),visc_sl(:,k));
%  end
% 
%  figure
% hold on
% plot(-convmfp(2:end)./Lx, zqm(2:end), '--b', 'LineWidth', 1)
% plot(-strmfp(2:end)./Lx,  zqm(2:end), ':b', 'LineWidth', 1)
% plot((-convmfp(2:end)-strmfp(2:end))./Lx, zqm(2:end), '-b', 'LineWidth', 1)
% plot((-viscmfp(2:end))./Lx, zqm(2:end), '-r', 'LineWidth', 1)
% plot((-strmfp(2:end)-convmfp(2:end)-viscmfp(2:end))./Lx, zqm(2:end), '-k', 'LineWidth', 1)
% % 