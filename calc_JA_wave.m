clear
close all
% baseDir = fullfile('/projects', 'standard', 'shenl','naras062','JA-Samvit','Retau-10','RUN-c0');
%baseDir = '/users/1/kuma0458/wave/wavy_wall';
%baseDir = '/users/1/kuma0458/wave/wave_c_2';
baseDir = '/users/1/kuma0458/wave/wavy_ret180';

%load('grid.mat')
tic
%ret=10;
ret=180;
nu=1/ret;
tstart=20200000;
step=  200000;
tend=  32000000;
% ff=fullfile(baseDir,'phi.mat');
fj=fullfile(baseDir,'flowrate.mat')
load(fj);
Ts = Jdot.*0;
phidots=Jdot.*0;
Tnls=Ts;
Tviscs=Ts;
check=Jdot.*0;
Nx=128;
Ny=16;
Nz=128;
%uphi=reshape(uphi,Nx,1,Nz);
%wphi=reshape(wphi,Nx,1,Nz);
%c=2;
c=0;
for tstep=tstart:step:tend
    % fn=		sprintf('DAT000%03d99999999',tstep)
    % fnmat= sprintf('gradflux000%03d99999999.mat',tstep)

    fn=sprintf('Sol%014d.h5',tstep);
    fnmat=sprintf('gradflux%014d.mat',tstep);
    fname = fullfile(baseDir,fn);
    fnamemat=fullfile(baseDir,fnmat);
    %fnmat = sprintf('grid%014d.mat',tstep);
    %fng = fullfile(baseDir,fnmat);
    fng=fullfile(baseDir,'grid.mat');
    load(fng);
    %fnp=sprintf('potvel%014d.mat',tstep);
fnp=sprintf('potvel.mat');
    fnp=fullfile(baseDir,fnp);
fnja = sprintf('jafields%014d.mat',tstep)
fnja = fullfile(baseDir,fnja);
load(fnp)
    % info=h5info(fname)
    fprintf('Reading %s\n', fname);
    u    = h5read(fname, '/u');
    v    = h5read(fname, '/v');
    w    = h5read(fname, '/w');
    time = h5read(fname, '/time')
    zz   = h5read(fname, '/zz');
    pey = h5read(fname,'/pey');
    pex = h5read(fname,'/pex');
    zz=zz';
    Ly=2*pi/pey;
    Lx=2*pi/pex;
    load(fnamemat)
    
    % xshift = c*time;
    % dx =X(3,1,1)-X(2,1,1);
    % ishift=round(xshift/dx);
    % uphis=circshift(uphi,ishift,1);
    % wphis=circshift(wphi,ishift,1);
uphis=uphi;
wphis=wphi;

    ox = dwdy-dvdz;
    oy = dudz-dwdx;
    oz = dvdx-dudy;
%[val id] = max(Z(:,1,1));
%xmax = X(id,1,1)
    % bot = 0.1*cos(2*(X(:,1,1)-c*time));
    % botshift = circshift(bot,-ishift,1);

         voz = single(v.*oz);
         woy = single(wc.*oy);
         uoy = single((u-c).*oy);
         vox = single(v.*ox);
         JAx = uphis.*(voz-woy+viscu);
         JAz = wphis.*(uoy-vox+viscw);
         JAnl = uphis.*(voz-woy)+wphis.*(uoy-vox);
         JAvisc = uphis.*viscu + wphis.*viscw;

clear u v w ox oy oz viscu viscw dwdy dvdz dudz dwdx dvdx dudy

        vol=JAx.*0+1;
        %%
        JAin=JAx+JAz;
JAnl(isnan(JAnl))=0;
JAvisc(isnan(JAvisc))=0;
JAtot = JAnl+JAvisc;
save(fnja,'JAnl','JAvisc');

    JAin(isnan(JAin))=0;
    JAin=squeeze(mean(JAin,2));
    JAin=trapz(zz',JAin,2);

    JAnlin = squeeze(mean(JAnl,2));
    JAviscin=squeeze(mean(JAvisc,2));
    JAnlin=trapz(zz',JAnlin,2);
    JAviscin=trapz(zz',JAviscin,2);


    Jacobian=1./dZetadz;
    T = -trapz(X(:,1,1),Jacobian.*JAin);
    Tnl = -trapz(X(:,1,1),Jacobian.*JAnlin);
    Tvisc = -trapz(X(:,1,1),Jacobian.*JAviscin);

   
    vol = trapz(X(:,1,1),Jacobian);
     it = find(t ==time)
   
 phidot = (Jdot(it))*Lx/Ly;
Ts(it)=T;
Tnls(it)=Tnl;
Tviscs(it)=Tvisc;
phidots(it)=phidot;


check(it)=100*((T+phidot)/vol);
end
%%
mja=fullfile(baseDir,'JAseries.mat')
save(mja,'Ts','Tnls','Tviscs','phidots','check')
%%
% save('JAseries.mat','Ts','Tnls','Tviscs','phidots','check')
cl=400;
%
figure
subplot(3,1,1)
pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(-JAnl(:,3,:)))
shading interp
c=colorbar
 clim([-cl cl])
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel(c,'$-u_{\phi}\cdot(u\times\omega )$','Interpreter','latex','FontSize',12)
ylabel('z/H')

subplot(3,1,2)
pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(-JAvisc(:,3,:)))
shading interp
c=colorbar
clim([-cl cl])
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
ylabel(c,'$-u_{\phi}\cdot(\nu \Delta u)$','Interpreter','latex','FontSize',12)
ylabel('z/H')

subplot(3,1,3)
pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(-JAvisc(:,3,:)-JAnl(:,3,:)))
shading interp
c=colorbar;
clim([-cl cl])
axis equal
xlim([0 2*pi])
ylim([-0.1 1])
xlabel('x/H')
ylabel('z/H')
ylabel(c,'$-u_{\phi}\cdot(u\times\omega +\nu \Delta u)$','Interpreter','latex','FontSize',12)
colormap jet
% %
% figure
% subplot(2,1,1)
% pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(uphis(:,1,:)))
% shading flat
% c=colorbar
% %clim([-2 2])
% axis equal
% xlim([0 2*pi])
% ylim([-0.1 1])
% ylabel(c,'$u_{\phi}$','Interpreter','latex','FontSize',12)
% ylabel('z/H')
% 
% subplot(2,1,2)
% pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(wphis(:,1,:)))
% shading flat
% c=colorbar
% %clim([-10 10])
% axis equal
% xlim([0 2*pi])
% xlabel('x/H')
% ylim([-0.1 1])
% ylabel(c,'$w_{\phi}$','Interpreter','latex','FontSize',12)
% ylabel('z/H')
% %%
% figure
% subplot(3,1,1)
% pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(u(:,1,:)))
% shading flat
% cc=colorbar;
% % clim([-1 1])
% axis equal
% xlim([0 2*pi])
% ylim([-0.1 1])
% ylabel(cc,'$u$','Interpreter','latex','FontSize',12)
% ylabel('z/H')
% 
% subplot(3,1,2)
% pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(w(:,1,:)))
% shading flat
% cc=colorbar;
% %clim([-10 10])
% axis equal
% xlim([0 2*pi])
% ylim([-0.1 1])
% ylabel(cc,'$w$','Interpreter','latex','FontSize',12)
% ylabel('z/H')
% 
% subplot(3,1,3)
% pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(uphi(:,1,:)))
% shading flat
% cc=colorbar;
% %clim([-10 10])
% axis equal
% xlim([0 2*pi])
% ylim([-0.1 1])
% ylabel(cc,'$\omega_y$','Interpreter','latex','FontSize',12)
% ylabel('z/H')
% xlabel('x/H')