clear
close all

baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
Nx=256;
Ny=192;
Nz=128;
uav=zeros(Nx,Ny,Nz);
vav=zeros(Nx,Ny,Nz);
wav=zeros(Nx,Ny,Nz);
oxav=zeros(Nx,Ny,Nz);
oyav=zeros(Nx,Ny,Nz);
ozav=zeros(Nx,Ny,Nz);
viscav=zeros(Nx,Ny,Nz);
nlav=zeros(Nx,Ny,Nz);
convav=zeros(Nx,Ny,Nz);
strav=zeros(Nx,Ny,Nz);
c=2;
pex=0.5;
kx=pex*[0:Nx/2-1,-Nx/2:-1]';
%kdis = exp(-(1i*ct).*kx);

%baseDir = '/users/1/kuma0458/wave/wavy_wall';
%baseDir = '/users/1/kuma0458/wave/wavy_ret180';
%baseDir = '/users/1/kuma0458/wave/wave_c_2';
%baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';
tstart=3025000000;
step  =   5000000;
tend  =3820000000;
counter=0;
for tstep=tstart:step:tend
	time = tstep/(1e+8)
	ct=c*time;
	kd = exp((1i*ct).*kx);
	kdis = reshape(kd,[Nx,1,1]);
	fn=sprintf('Sol%014d.h5',tstep);
	fnmat=sprintf('gradflux%014d.mat',tstep);
	fname = fullfile(baseDir,fn);
	u    = h5read(fname, '/u')-c;
	v    = h5read(fname, '/v');

	fnamemat=fullfile(baseDir,fnmat);
        mt=matfile(fnamemat);
	w=mt.wc;
	ox = mt.dwdy-mt.dvdz;
	oy = mt.dudz-mt.dwdx;
	oz = mt.dvdx-mt.dudy;
	clear mt;
	
	fngs = sprintf('jafields%014d.mat',tstep)
 	fng = fullfile(baseDir,fngs);
 	load(fng)
        
	fvisc = fft(JAvisc,[],1).*kdis;
        viscav = viscav+ifft(fvisc,[],1,'symmetric');
        clear fvisc JAvisc
        fnl = fft(JAnl,[],1).*kdis;
        nlav = nlav+ifft(fnl,[],1,'symmetric');
        clear fnl JAnl	
	fconv = fft(JAconv,[],1).*kdis;
        convav = convav+ifft(fconv,[],1,'symmetric');
        clear fconv JAconv
	fstr = fft(JAstr,[],1).*kdis;
        strav = strav+ifft(fstr,[],1,'symmetric');
        clear fstr JAstr

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
viscav=viscav./counter;
nlav = nlav./counter;
convav = convav./counter;
strav = strav./counter;

mf=fullfile(baseDir,'mean_fields.mat');
save(mf,'uav','vav','wav','oxav','oyav','ozav','viscav','nlav','strav','convav');
