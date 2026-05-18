clear
close all
baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';

Nx=128;
Ny=32;
Nz=128;
uav=zeros(Nx,Ny,Nz);
vav=zeros(Nx,Ny,Nz);
wav=zeros(Nx,Ny,Nz);
oxav=zeros(Nx,Ny,Nz);
oyav=zeros(Nx,Ny,Nz);
ozav=zeros(Nx,Ny,Nz);
c=2;
pex=1;
kx=pex*[0:Nx/2-1,0,-Nx/2+1:-1]';
%kdis = exp(-(1i*ct).*kx);

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
	time = tstep/(1e+8)
	ct=c*time;
	kd = exp((1i*ct).*kx);
	kdis = reshape(kd,[Nx,1,1]);
	fn=sprintf('Sol%014d.h5',tstep);
	fnmat=sprintf('gradflux%014d.mat',tstep);
	fname = fullfile(baseDir,fn);
	fnamemat=fullfile(baseDir,fnmat);
	mt=matfile(fnamemat);
	u    = h5read(fname, '/u')-c;
	v    = h5read(fname, '/v');
	w=mt.wc;
	ox = mt.dwdy-mt.dvdz;
	oy = mt.dudz-mt.dwdx;
	oz = mt.dvdx-mt.dudy;
	clear mt;
	
	uav = uav + ifft( (fft(u,[],1).*kdis),[],1,'symmetric');
	vav = vav + ifft( (fft(v,[],1).*kdis),[],1,'symmetric');
	wav = wav + ifft( (fft(w,[],1).*kdis),[],1,'symmetric');
	oxav= oxav+ ifft( (fft(ox,[],1).*kdis),[],1,'symmetric');
	oyav= oyav+ ifft( (fft(oy,[],1).*kdis),[],1,'symmetric');
	ozav= ozav+ ifft( (fft(oz,[],1).*kdis),[],1,'symmetric');
	counter =counter+1;
end
uav=uav./counter;
vav=vav./counter;
wav=wav./counter;
oxav=oxav./counter;
oyav=oyav./counter;
ozav=ozav./counter;

save('mean_fields.mat','uav','vav','wav','oxav','oyav','ozav');
