clear
close all
%baseDir = '/users/1/kuma0458/wave/wavy_ret180';
%baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';
%baseDir = '/users/1/kuma0458/wave/c2ak2_re180/run';
baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
Nx=256;
Ny=192;
Nz=128;
wave_n=12;
fn    = 'grid.h5';
fname = fullfile(baseDir,fn);
pex  = h5read(fname, '/pex');
kx=pex*[0:Nx/2-1,-Nx/2:-1]';
Lx=2*pi/pex;
x=[0:Nx-1]*Lx/Nx;

upav=zeros(Nx,Ny,Nz);
vpav=zeros(Nx,Ny,Nz);
wpav=zeros(Nx,Ny,Nz);
oxpav=zeros(Nx,Ny,Nz);
oypav=zeros(Nx,Ny,Nz);
ozpav=zeros(Nx,Ny,Nz);
viscpav=zeros(Nx,Ny,Nz);
nlpav=zeros(Nx,Ny,Nz);
convpav=zeros(Nx,Ny,Nz);
strpav=zeros(Nx,Ny,Nz);

mf=fullfile(baseDir,'mean_fields.mat');
 load(mf)
for id = 0:wave_n-1
	xdis  = id*L0;
	kd    = exp((1i*xdis).*kx);
	kdis  = reshape(kd,[Nx,1,1]);
	upav  = upav  + ifft( (fft(uav,[],1).*kdis),[],1,'symmetric');
	vpav  = vpav  + ifft( (fft(vav,[],1).*kdis),[],1,'symmetric');
	wpav  = wpav  + ifft( (fft(wav,[],1).*kdis),[],1,'symmetric');
	oxpav = oxpav + ifft( (fft(oxav,[],1).*kdis),[],1,'symmetric');
	oypav = oypav + ifft( (fft(oyav,[],1).*kdis),[],1,'symmetric');
	ozpav = ozpav + ifft( (fft(ozav,[],1).*kdis),[],1,'symmetric');
	viscpav = viscpav + ifft((fft(viscav,[],1).*kdis),[],1,'symmetric');
	  nlpav =   nlpav + ifft((fft(  nlav,[],1).*kdis),[],1,'symmetric');
	  convpav = convpav + ifft((fft(convav,[],1).*kdis),[],1,'symmetric');
          strpav =   strpav + ifft((fft(  strav,[],1).*kdis),[],1,'symmetric');
end

upav = squeeze(mean(upav(1:Nxl,:,:),2))./wave_n;
vpav = squeeze(mean(vpav(1:Nxl,:,:),2))./wave_n;
wpav = squeeze(mean(wpav(1:Nxl,:,:),2))./wave_n;

oxpav =squeeze(mean(oxpav(1:Nxl,:,:),2))./wave_n;
oypav =squeeze(mean(oypav(1:Nxl,:,:),2))./wave_n;
ozpav =squeeze(mean(ozpav(1:Nxl,:,:),2))./wave_n;

viscpav =squeeze(mean(viscpav(1:Nxl,:,:),2))./wave_n;
nlpav =  squeeze(mean(  nlpav(1:Nxl,:,:),2))./wave_n;
convpav =  squeeze(mean(convpav(1:Nxl,:,:),2))./wave_n;
strpav =  squeeze(mean( strpav(1:Nxl,:,:),2))./wave_n;

mf=fullfile(baseDir,'mean_phase_fields.mat');
save(mf,'upav','vpav','wpav','oxpav','oypav','ozpav','viscpav','nlpav','convpav','strpav');

