clear
close all
%baseDir = '/users/1/kuma0458/wave/wavy_ret180';
%baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';
%baseDir = '/users/1/kuma0458/wave/c2ak2_re180/run';
baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
Nx=256;
Ny=192;
Nz=128;
ak=0.2;
c=2;
convp=zeros(Nx,Ny,Nz);
convn=zeros(Nx,Ny,Nz);
strp=zeros(Nx,Ny,Nz);
strn=zeros(Nx,Ny,Nz);

wave_n=12;
fname = fullfile(baseDir,'grid.h5');
z   = h5read(fname, '/zz');
zw   = h5read(fname, '/zw');
pex  = h5read(fname, '/pex');
pey = h5read(fname, '/pey');
Ly=2*pi/pey;
Lx=2*pi/pex;
L0=Lx/wave_n;
kx=pex*[0:Nx/2-1,-Nx/2:-1]';
Nxl = round(Nx/wave_n)+1;
load(fullfile(baseDir,'phi_interp_2d.mat'),'uphi','wphi')
up=sqrt(uphi.^2+wphi.^2);
uph=uphi./up;
wph=wphi./up;
% ai = uph;  ak = wph;
ci = -wph; ck = uph;
clear uph uphi wph wphi up

tstart=3025000000;
step  =   5000000;
tend  =3820000000;
counter=0;
for tstep=tstart:step:tend

counter=counter+1;

    fn=sprintf('Sol%014d.h5',tstep);
    fname = fullfile(baseDir,fn);
    time=h5read(fname,'/time');
    %time = tstep/(1e+8)
    ct=c*time;
    kd = exp((1i*ct).*kx);
    kdis = reshape(kd,[Nx,1,1]);
    u= h5read(fname, '/u')-c;
    u=ifft( (fft(u,[],1).*kdis),[],1,'symmetric');

    fnmat=sprintf('gradflux%014d.mat',tstep);
    fnamemat=fullfile(baseDir,fnmat);
    mt=matfile(fnamemat);
    w=mt.wc;
    w=ifft( (fft(w,[],1).*kdis),[],1,'symmetric');

    uc = ci.*u + ck.*w;
    mp = double(uc>0);
    mn = double(uc<0);

    fngs = sprintf('jafields%014d.mat',tstep)
        fng = fullfile(baseDir,fngs);
        load(fng)
	
	%fconv = fft(JAconv,[],1).*kdis;
        %conv = ifft(fconv,[],1,'symmetric');	
       	convp=convp+mp.*JAconv;
	convn=convn+mn.*JAconv;

	%fstr = fft(JAstr,[],1).*kdis;
        %str = ifft(fstr,[],1,'symmetric');
         strp = strp + mp.*JAstr;
    	 strn = strn + mn.*JAstr;

end
convp =squeeze(mean(convp./counter,2));
convn =squeeze(mean(convn./counter,2));
strp = squeeze(mean( strp./counter,2));
strn = squeeze(mean( strn./counter,2));

%for id = 0:wave_n-1
%        xdis  = id*L0;
%        kd    = exp((1i*xdis).*kx);
%        kdis  = reshape(kd,[Nx,1,1]);
%        convp  = convp  + ifft( (fft(convp,[],1).*kdis),[],1,'symmetric');
%        convn  = convn  + ifft( (fft(convn,[],1).*kdis),[],1,'symmetric');
%        strp  = strp  + ifft( (fft(strp,[],1).*kdis),[],1,'symmetric');
%        strn  = strn  + ifft( (fft(strn,[],1).*kdis),[],1,'symmetric');
%end
%convp = squeeze(mean(convp(1:Nxl,:,:),2))./wave_n;
%convn = squeeze(mean(convn(1:Nxl,:,:),2))./wave_n;
%strp = squeeze(mean(strp(1:Nxl,:,:),2))./wave_n;
%strn = squeeze(mean(strn(1:Nxl,:,:),2))./wave_n;


fp = fullfile(baseDir,'partial_sum.mat');
save(fp,'convn','convp','strn','strp')

