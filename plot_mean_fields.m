close all
clear
baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
Nx=256;
Ny=192;
Nz=128;
wave_n=12;
nc=15;
load(fullfile( baseDir,'grid.mat'))
load(fullfile( baseDir,'partial_sum.mat'))
Nxl = round(Nx/wave_n)+1;
xc = squeeze(X(1:Nxl,1,:));
zc = squeeze(Z(1:Nxl,1,:));
cl=5e+3;
figure
subplot(1,4,1)
contourf(xc,zc,-convp,nc)
axis equal
xlim([0 xc(end,end)])
ylim([-0.04 1])
colorbar 
clim([-cl cl])

subplot(1,4,2)
contourf(xc,zc,-strp,nc)
axis equal
xlim([0 xc(end,end)])
ylim([-0.04 1])
colorbar 
clim([-cl cl])

subplot(1,4,3)
contourf(xc,zc,-convn,nc)
axis equal
xlim([0 xc(end,end)])
ylim([-0.04 1])
colorbar 
clim([-cl cl])

subplot(1,4,4)
contourf(xc,zc,-strn,nc)
axis equal
xlim([0 xc(end,end)])
ylim([-0.04 1])
colorbar 
clim([-cl cl])

colormap jet
%%
cl1=80;
load(fullfile(baseDir,'mean_phase_fields.mat'))
figure
subplot(1,3,1)
contourf(xc,zc,-nlpav,nc)
axis equal
xlim([0 xc(end,end)])
ylim([-0.04 1])
clim([-cl1 cl1])
colorbar 

subplot(1,3,2)
contourf(xc,zc,-viscpav,nc)
axis equal
xlim([0 xc(end,end)])
ylim([-0.04 1])
clim([-cl1 cl1])
colorbar 

subplot(1,3,3)
contourf(xc,zc,-viscpav-nlpav,nc)
axis equal
xlim([0 xc(end,end)])
ylim([-0.04 1])
clim([-cl1 cl1])
colorbar 
colormap jet