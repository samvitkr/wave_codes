%% Read timestep 10 from OCF data and save into a .mat file
clear
close all
%baseDir = '/users/1/kuma0458/wave/wave_c_2';
%baseDir = '/users/1/kuma0458/wave/wave_c_2';
%baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';

baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
c=2;

%baseDir = '/scratch.global/kuma0458/c14ak1_re180/run';
%c=14;

 % baseDir = '/scratch.global/kuma0458/c8ak1_re180/run';
 % c=8;

Nx=256;
Ny=192;
Nz=128;

x1=10;
y1=10;
width=1098;
height=541;

% tstart= 4300000000;
%   step=    1250000;
%  tend = 4770000000;

tstart=3025000000;
step  =   5000000;
tend  =4020000000;

fnmat='grid.mat';
fng = fullfile(baseDir,fnmat);
load(fng)

load(fullfile(baseDir,'phi_interp_2d.mat'));%,'uphi','wphi','J')



 yl=-0.05;
ret=180;
cnl=50;
cv=50;
snum=15;
%%
ts=[tstart:step:tend];
tstep=ts(snum);
fnmat=sprintf('jafields%014d.mat',tstep);
    fnamemat = fullfile(baseDir,fnmat);
    load(fnamemat);
   fn=sprintf('Sol%014d.h5',tstep)
fname   = fullfile(baseDir,fn);
time = h5read(fname ,'/time')
sgt=sprintf("t=%0.4f",time);




%%
figure
t=tiledlayout(2,1);
nexttile
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(-mean(JAstr(:,:,:),2)));
shading flat
axis equal
xlim([0 2*pi])
ylim([yl 1])
ylabel('z/H')
cu=colorbar;
ylabel(cu,'stretching tilting','Interpreter','latex','FontSize',12)
 %clim([-3 3])
clim([-cnl cnl])

nexttile
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(-mean(JAconv(:,:,:),2)));
shading flat
axis equal
xlim([0 2*pi])
ylim([yl 1])
ylabel('z/H')
%clim([-30 30])
clim([-cv cv])
cw=colorbar;
ylabel(cw,'convective','Interpreter','latex','FontSize',12)


% nexttile
% pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(-mean(JAnl(:,:,:),2 )));
% shading flat
% axis equal
% xlim([0 2*pi])
% ylim([yl 1])
% ylabel('z/H')
% co=colorbar;
% clim([-cnl cnl])
% %clim([-30 30])
% ylabel(co,'$-u_{\phi}\cdot(u\times\omega+\nu\Delta u)$','Interpreter','latex','FontSize',12)
% xlabel('x/H')
sgtitle(sgt)
colormap jet


%%
figure
t=tiledlayout(3,1);
nexttile
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(-mean(JAnl(:,:,:),2)));
shading flat
axis equal
xlim([0 2*pi])
ylim([yl 1])
ylabel('z/H')
cu=colorbar;
ylabel(cu,'$-u_{\phi}\cdot(u\times\omega)$','Interpreter','latex','FontSize',12)
 %clim([-3 3])
clim([-cnl cnl])

nexttile
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(-mean(JAvisc(:,:,:),2)));
shading flat
axis equal
xlim([0 2*pi])
ylim([yl 1])
ylabel('z/H')
%clim([-30 30])
clim([-cv cv])
cw=colorbar;
ylabel(cw,'$-u_{\phi}\cdot(\nu \Delta u)$','Interpreter','latex','FontSize',12)

JA=JAvisc+JAnl;
nexttile
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(-mean(JA(:,:,:),2 )));
shading flat
axis equal
xlim([0 2*pi])
ylim([yl 1])
ylabel('z/H')
co=colorbar;
clim([-cnl cnl])
%clim([-30 30])
ylabel(co,'$-u_{\phi}\cdot(u\times\omega+\nu\Delta u)$','Interpreter','latex','FontSize',12)
xlabel('x/H')
sgtitle(sgt)
colormap jet
%%
load(fullfile(baseDir,'JAseries_statwave.mat'))
figure
plot(check,'-k')
xline(snum)
%%
cw=max(abs(wphi),[],'all');

figure
subplot(3,1,1)
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze( wphi(:,1,:)));
shading flat
colorbar
subplot(3,1,2)
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze( uphi(:,1,:)));
shading flat
colorbar

subplot(3,1,3)
contourf(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),( Phi_2D(:,:)),120);
shading flat
colorbar
