clear
close all
% baseDir = fullfile('/projects', 'standard', 'shenl','naras062','JA-Samvit','Retau-10','RUN-c0');
baseDir = '/users/1/kuma0458/wave/wavy_wall';
load('grid.mat')
tic
ret=10;
nu=1/ret;
tstart=2000000;
step=2000000;
tend=2000000;
ff=fullfile(baseDir,'pot.mat');
fj=fullfile(baseDir,'flowrate.mat')
load(ff)
load(fj);
Nx=128;
Ny=16;
Nz=128;
uphi=reshape(uphi,Nx,1,Nz);
wphi=reshape(wphi,Nx,1,Nz);
load()

for tstep=tstart:step:tend
    % fn=		sprintf('DAT000%03d99999999',tstep)
    % fnmat= sprintf('gradflux000%03d99999999.mat',tstep)

    fn=sprintf('Sol%014d.h5',tstep);
    fnmat=sprintf('gradflux%014d.mat',tstep);
    fname = fullfile(baseDir,fn);
    fnamemat=fullfile(baseDir,fnmat);
    % info=h5info(fname)
    fprintf('Reading %s\n', fname);
    u    = h5read(fname, '/u');
    v    = h5read(fname, '/v');
        w    = h5read(fname, '/w');

    zz   = h5read(fname, '/zz');
    zz=zz';
    load(fnamemat)
    ox = dwdy-dvdz;
    oy = dudz-dwdx;
    oz = dvdx-dudy;

    voz = single(v.*oz);
    woy = single(wc.*oy);
    uoy = single(u.*oy);
    vox = single(v.*ox);
    JAx = uphi.*(voz-woy+viscu);
    JAz = wphi.*(uoy-vox+viscw);
    JAnl = uphi.*(voz-woy)+wphi.*(uoy-vox);
    JAvisc = uphi.*viscu + wphi.*viscw;

    vol=JAx.*0+1;
    %%
    JAin=JAx+JAz;
    
JAin(isnan(JAin))=0;
JAin=squeeze(mean(JAin,2));
JAin=trapz(zz',JAin,2);

Jacobian=1./dZetadz;
T = -trapz(X(:,1,1),Jacobian.*JAin);
vol = trapz(X(:,1,1),Jacobian);

end
%%
figure
subplot(3,1,1)
pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(-JAnl(:,3,:)))
shading flat
c=colorbar
clim([-2 2])
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel(c,'$-u_{\phi}\cdot(u\times\omega )$','Interpreter','latex','FontSize',12)
ylabel('z/H')

subplot(3,1,2)
pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(-JAvisc(:,3,:)))
shading flat
c=colorbar
clim([-10 10])
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel(c,'$-u_{\phi}\cdot(\nu \Delta u)$','Interpreter','latex','FontSize',12)
ylabel('z/H')

subplot(3,1,3)
pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(-JAvisc(:,3,:)-JAnl(:,3,:)))
shading flat
c=colorbar;
clim([-10 10])
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
xlabel('x/H')
ylabel('z/H')
ylabel(c,'$-u_{\phi}\cdot(u\times\omega +\nu \Delta u)$','Interpreter','latex','FontSize',12)
colormap jet
%%
figure
subplot(2,1,1)
pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(uphi(:,1,:)))
shading flat
c=colorbar
%clim([-2 2])
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel(c,'$u_{\phi}$','Interpreter','latex','FontSize',12)
ylabel('z/H')

subplot(2,1,2)
pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(wphi(:,1,:)))
shading flat
c=colorbar
%clim([-10 10])
axis equal
xlim([0 2*pi])
xlabel('x/H')
ylim([-0.1 1])
ylabel(c,'$w_{\phi}$','Interpreter','latex','FontSize',12)
ylabel('z/H')
%%
figure
subplot(3,1,1)
pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(u(:,1,:)))
shading flat
c=colorbar
%clim([-2 2])
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel(c,'$u$','Interpreter','latex','FontSize',12)
ylabel('z/H')

subplot(3,1,2)
pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(w(:,1,:)))
shading flat
c=colorbar
%clim([-10 10])
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel(c,'$w$','Interpreter','latex','FontSize',12)
ylabel('z/H')

subplot(3,1,3)
pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(oy(:,1,:)))
shading flat
c=colorbar
%clim([-10 10])
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel(c,'$\omega_y$','Interpreter','latex','FontSize',12)
ylabel('z/H')
xlabel('x/H')