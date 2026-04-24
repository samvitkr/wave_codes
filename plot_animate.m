%% Read timestep 10 from OCF data and save into a .mat file
clear
close all
%baseDir = '/users/1/kuma0458/wave/wave_c_2';
%baseDir = '/users/1/kuma0458/wave/wave_c_2';
baseDir = '/users/1/kuma0458/wave/wavy_ret180';

x1=10;
y1=10;
width=1096;
height=582;
tstart=20200000;
step=200000;
tend=32000000;
ret=180;
%%
fvn=fullfile(baseDir,'vel_field.avi');
%v=VideoWriter('vel_field',);
v = VideoWriter(fvn, 'Motion JPEG AVI');
v.FrameRate=2;
open(v)
f=figure('OuterPosition',[x1 y1 width height]);
for tstep=tstart:step:tend
% fn=sprintf('DAT000%03d99999999',tstep)
fn=sprintf('Sol%014d.h5',tstep)
fname   = fullfile(baseDir,fn);   % <-- edit if naming differs
% info=h5info(fname)
%fnmat = sprintf('grid%014d.mat',tstep);
fnmat='grid.mat';
fng = fullfile(baseDir,fnmat);

fnmat=sprintf('gradflux%014d.mat',tstep);
    fnamemat = fullfile(baseDir,fnmat);
    m=matfile(fnamemat);
        oy = m.dudz-m.dwdx;

u    = h5read(fname, '/u');
w    = h5read(fname, '/w');
time = h5read(fname ,'/time')
sgt=sprintf("t=%0.4f",time);

load(fng)
%%

t=tiledlayout(3,1);
nexttile
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(u(:,1,:)));
shading flat
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel('z/H')
cu=colorbar;
ylabel(cu,'u^+','FontSize',12)
%clim([-2 10])
clim([-5 20])

nexttile
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(w(:,1,:)));
shading flat
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel('z/H')
%clim([-0.5 0.5])
clim([-4 4])
cw=colorbar;
ylabel(cw,'w^+','FontSize',12)

nexttile
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(oy(:,1,:))./ret);
shading flat
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel('z/H')
co=colorbar;
%clim([-0.5 2.5])
clim([-4 4])
ylabel(co,'\omega_y^+','FontSize',12)
xlabel('x/H')
sgtitle(sgt)
frame=getframe(f);
writeVideo(v,frame);
clf
end

%%
cnl=400;
cv=200;
close all
fvja=fullfile(baseDir,"Ja_field.avi")
 %vv = VideoWriter('Ja_field.avi', 'Motion JPEG AVI');
  vv = VideoWriter(fvja, 'Motion JPEG AVI');

 vv.FrameRate=2;
 open(vv)
f=figure('OuterPosition',[x1 y1 width height]);
for tstep=tstart:step:tend
% fn=sprintf('DAT000%03d99999999',tstep)
fn=sprintf('Sol%014d.h5',tstep)
fname   = fullfile(baseDir,fn);   % <-- edit if naming differs
% info=h5info(fname)
%fnmat = sprintf('grid%014d.mat',tstep);
fnmat='grid.mat';
fng = fullfile(baseDir,fnmat);

fnmat=sprintf('jafields%014d.mat',tstep);
    fnamemat = fullfile(baseDir,fnmat);
    load(fnamemat);
   
time = h5read(fname ,'/time')
sgt=sprintf("t=%0.4f",time);

load(fng)
%

t=tiledlayout(3,1);
nexttile
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(-JAnl(:,1,:)));
shading flat
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel('z/H')
cu=colorbar;
ylabel(cu,'$-u_{\phi}\cdot(u\times\omega)$','Interpreter','latex','FontSize',12)
 %clim([-3 3])
clim([-cnl cnl])

nexttile
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(-JAvisc(:,1,:)));
shading flat
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel('z/H')
%clim([-30 30])
clim([-cv cv])
cw=colorbar;
ylabel(cw,'$-u_{\phi}\cdot(\nu \Delta u)$','Interpreter','latex','FontSize',12)

JA=JAvisc+JAnl;
nexttile
pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(-JA(:,1,:)));
shading flat
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel('z/H')
co=colorbar;
clim([-cnl cnl])
%clim([-30 30])
ylabel(co,'$-u_{\phi}\cdot(u\times\omega+\nu\Delta u)$','Interpreter','latex','FontSize',12)
xlabel('x/H')
sgtitle(sgt)
colormap jet
 frame=getframe(f);
 writeVideo(vv,frame);
 clf
end
%%
% close all
 vol=2*pi;
fla=fullfile(baseDir,'JAseries.mat')
load(fla)
% load("JAseries.mat")
fj=fullfile(baseDir,'flowrate.mat')
load(fj)
% load("../wave_c_2/flowrate.mat")
subplot(1,2,1)
hold on
plot(t,Tnls./vol,'--b','LineWidth',1.5)
plot(t,Tviscs./vol,':b','LineWidth',1.5)
plot(t,Ts./vol,'-b','LineWidth',1.5)
plot(t,phidots./vol,'-r','LineWidth',1.5)
plot(t,(Ts+phidots)./vol,'-k','LineWidth',1.5)
yline(1)

hold off
legend('T_{NL}','T_{visc}','T','\phi_t','T+\phi_t','Location','eastoutside')
legend boxoff
xlabel('t')
subplot(1,2,2)

plot(t,100-check,'-.k')
xlabel('t')
ylabel('%Error')
