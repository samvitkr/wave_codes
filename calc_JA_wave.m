clear
close all
% baseDir = fullfile('/projects', 'standard', 'shenl','naras062','JA-Samvit','Retau-10','RUN-c0');
%baseDir = '/users/1/kuma0458/wave/wavy_wall';
%baseDir = '/users/1/kuma0458/wave/wave_c_2';
%baseDir = '/users/1/kuma0458/wave/wavy_ret180';
%baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';
tic
%ret=10;
ret=180;
nu=1/ret;
baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
Nx=256;
Ny=192;
Nz=128;
fname = fullfile(baseDir,'grid.h5');
zz   = h5read(fname, '/zz');
zw   = h5read(fname, '/zw');
pex  = h5read(fname, '/pex');
pey = h5read(fname, '/pey');
Ly=2*pi/pey;
Lx=2*pi/pex;
kx=pex*[0:Nx/2-1,-Nx/2:-1]';
load(fullfile(baseDir,'flowrate.mat'))
load(fullfile(baseDir,'grid.mat'))
load(fullfile(baseDir,'phi_interp_2d.mat'),'uphi','wphi')

fup=fft(uphi,[],1);
fwp=fft(wphi,[],1);
clear uphi wphi

Ts = Jdot.*0;
phidots=Jdot.*0;
Tnls=Ts;
Tviscs=Ts;
Tconvs=Ts;
Tstrs=Ts;
check=Jdot.*0;
c=2;

tstart=3025000000;
step  =   5000000;
tend  =3820000000;

for tstep=tstart:step:tend
	fn=sprintf('Sol%014d.h5',tstep);
	fname = fullfile(baseDir,fn);
	fnmat=sprintf('gradflux%014d.mat',tstep);
	load(fullfile(baseDir,fnmat));
	fngn = sprintf('grid%014d.mat',tstep);
	load(fullfile(baseDir,fngn));
	
	%fnp=sprintf('potvel%014d.mat',tstep);
	%fnp=sprintf('potvel.mat');
	%fnp=fullfile(baseDir,fnp);
	fnja = sprintf('jafields%014d.mat',tstep)
	fnja = fullfile(baseDir,fnja);
	%load(fnp)
	% info=h5info(fname)
	fprintf('Reading %s\n', fname);
	u    = h5read(fname, '/u');
	v    = h5read(fname, '/v');
	time = h5read(fname, '/time')
	ct=c*time;
	
	kd = exp((-1i*ct).*kx);
	kdis = reshape(kd,[Nx,1,1]);
	uphis=ifft( (fup.*kdis),[],1,'symmetric');
	wphis=ifft( (fwp.*kdis),[],1,'symmetric');
	
	ox = dwdy-dvdz;
	oy = dudz-dwdx;
	oz = dvdx-dudy;
	
	voz = single(v.*oz);
	woy = single(wc.*oy);
	uoy = single((u-c).*oy);
	vox = single(v.*ox);
	
    % JAx = uphis.*(voz-woy+viscu);
	% JAz = wphis.*(uoy-vox+viscw);

	JAnl = uphis.*(voz-woy)+wphis.*(uoy-vox);
	JAconv=uphis.*(   -woy)+wphis.*(uoy    );
    JAstr =uphis.*(voz    )+wphis.*(   -vox);

    JAvisc = uphis.*viscu + wphis.*viscw;


	clear u v w ox oy oz viscu viscw dwdy dvdz dudz dwdx dvdx dudy
	
	vol=JAnl.*0+1;
	%%
	% JAin=JAx+JAz;
	JAnl(isnan(JAnl))=0;
	JAvisc(isnan(JAvisc))=0;
    JAconv(isnan(JAconv))=0;
    JAstr(isnan(JAstr))=0;

	JAtot = JAnl+JAvisc;

	save(fnja,'JAnl','JAvisc','JAconv','JAstr');
	
    JAin=JAtot;
	JAin(isnan(JAin))=0;
	JAin=squeeze(mean(JAin,2));
	JAin=trapz(zz,JAin,2);
	
	JAnlin = squeeze(mean(JAnl,2));
	JAviscin=squeeze(mean(JAvisc,2));
	JAnlin=trapz(zz,JAnlin,2);
	JAviscin=trapz(zz,JAviscin,2);
	
    JAconv = squeeze(mean(JAconv,2));
	JAconv = trapz(zz,JAconv,2);
	
    JAstr = squeeze(mean(JAstr,2));
	JAstr = trapz(zz,JAstr,2);

	Jacobian=1./dZetadz;

	T = -trapz(X(:,1,1),Jacobian.*JAin);
	Tnl = -trapz(X(:,1,1),Jacobian.*JAnlin);
	Tvisc = -trapz(X(:,1,1),Jacobian.*JAviscin);
	Tconv= -trapz(X(:,1,1),Jacobian.*JAconv);
    Tstr = -trapz(X(:,1,1),Jacobian.*JAstr);

	vol = trapz(X(:,1,1),Jacobian);
	it = find(t ==time)
	phidot = (Jdot(it))*Lx/Ly;
	Ts(it)=T;
	  Tnls(it)=Tnl;
	Tviscs(it)=Tvisc;
    Tconvs(it)=Tconv;
     Tstrs(it)=Tstr;
	phidots(it)=phidot;
	check(it)=100*((T+phidot)/vol);
end
%%
mja=fullfile(baseDir,'JAseries.mat')
save(mja,'Ts','Tnls','Tviscs','Tconvs','Tstrs','phidots','check')
%%
% cl=400;
% %
% figure
% subplot(3,1,1)
% pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(-JAnl(:,3,:)))
% shading interp
% c=colorbar
% clim([-cl cl])
% axis equal
% xlim([0 2*pi])
% ylim([-0.1 1])
% ylabel(c,'$-u_{\phi}\cdot(u\times\omega )$','Interpreter','latex','FontSize',12)
% ylabel('z/H')
% 
% subplot(3,1,2)
% pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(-JAvisc(:,3,:)))
% shading interp
% c=colorbar
% clim([-cl cl])
% axis equal
% xlim([0 2*pi])
% ylim([-0.1 1])
% ylabel(c,'$-u_{\phi}\cdot(\nu \Delta u)$','Interpreter','latex','FontSize',12)
% ylabel('z/H')
% 
% subplot(3,1,3)
% pcolor(squeeze(X(:,3,:)),squeeze(Z(:,3,:)),squeeze(-JAvisc(:,3,:)-JAnl(:,3,:)))
% shading interp
% c=colorbar;
% clim([-cl cl])
% axis equal
% xlim([0 2*pi])
% ylim([-0.1 1])
% xlabel('x/H')
% ylabel('z/H')
% ylabel(c,'$-u_{\phi}\cdot(u\times\omega +\nu \Delta u)$','Interpreter','latex','FontSize',12)
% colormap jet
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
