clear
close all
tic
%ret=10;
ret=180;
nu=1/ret;
%baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
%c=2;
baseDir = '/scratch.global/kuma0458/c14ak1_re180/run';
c=14;
%baseDir = '/scratch.global/kuma0458/c24ak2_re180/run';
%c=24;

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
load(fullfile(baseDir,'phi_interp_2d.mat'),'uphi','wphi','J')
	Jacobian=1./dZetadz;
	vol = trapz(X(:,1,1),Jacobian);
uphi=uphi./J;
wphi = wphi./J;
% fup=fft(uphi,[],1);
% fwp=fft(wphi,[],1);
% clear uphi wphi

Ts = Jdot.*0;
phidots=Jdot.*0;
Tnls=Ts;
Tviscs=Ts;
Tconvs=Ts;
Tstrs=Ts;
check=Jdot.*0;

%tstart=3025000000;
%step  =   5000000;
%tend  =3820000000;

% tstart=3020400000;
% step  =    400000;
% tend  =3084000000;

% tstart=3020200000;
% step  =    200000;
% tend  =3052000000;

tstart=3550000000;
step  =    400000;
tend  =3578000000;

for tstep=tstart:step:tend
	fn=sprintf('Sol%014d.h5',tstep);
	fname = fullfile(baseDir,fn);
	fnmat=sprintf('gradflux%014d.mat',tstep);
	load(fullfile(baseDir,fnmat));
	%fngn = sprintf('grid%014d.mat',tstep);
	%load(fullfile(baseDir,fngn));
	

	 fnja = sprintf('jafields%014d.mat',tstep)
	 fnja = fullfile(baseDir,fnja);

	fprintf('Reading %s\n', fname);
	u    = h5read(fname, '/u');
	v    = h5read(fname, '/v');
	time = h5read(fname, '/time')
	ct=c*time;
	
	% kd = exp((-1i*ct).*kx);
    kd = exp((1i*ct).*kx);
	kdis = reshape(kd,[Nx,1,1]);

    u = ifft( (fft(u,[],1).*kdis),[],1,'symmetric');
	v = ifft( (fft(v,[],1).*kdis),[],1,'symmetric');
	wc= ifft( (fft(wc,[],1).*kdis),[],1,'symmetric');	
	ox = dwdy-dvdz;
	oy = dudz-dwdx;
	oz = dvdx-dudy;
    clear dwdy dvdz dudz dwdx dvdx dudy
	ox= ifft( (fft(ox,[],1).*kdis),[],1,'symmetric');
	oy= ifft( (fft(oy,[],1).*kdis),[],1,'symmetric');
	oz= ifft( (fft(oz,[],1).*kdis),[],1,'symmetric');
    
	voz = single(v.*oz);
	woy = single(wc.*oy);
	uoy = single((u-c).*oy);
	vox = single(v.*ox);
	clear u v wc ox oy oz
    viscu= ifft( (fft(viscu,[],1).*kdis),[],1,'symmetric');
    viscw= ifft( (fft(viscw,[],1).*kdis),[],1,'symmetric');

	JAnl = uphi.*(voz-woy)+wphi.*(uoy-vox);
	JAconv=uphi.*(   -woy)+wphi.*(uoy    );
    JAstr =uphi.*(voz    )+wphi.*(   -vox);

    JAvisc = uphi.*viscu + wphi.*viscw;


	clear u v w ox oy oz viscu viscw dwdy dvdz dudz dwdx dvdx dudy
	
	%%
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


	T = -trapz(X(:,1,1),Jacobian.*JAin);
	Tnl = -trapz(X(:,1,1),Jacobian.*JAnlin);
	Tvisc = -trapz(X(:,1,1),Jacobian.*JAviscin);
	Tconv= -trapz(X(:,1,1),Jacobian.*JAconv);
    Tstr = -trapz(X(:,1,1),Jacobian.*JAstr);

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
mja=fullfile(baseDir,'JAseries_statwave.mat')
save(mja,'Ts','Tnls','Tviscs','Tconvs','Tstrs','phidots','check')
