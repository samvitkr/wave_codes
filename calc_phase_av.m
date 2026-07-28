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
lambda = 2*pi./(wave_n*pex);

for ii=1:wave_n-1
    if(mod(ii*Nx,wave_n)==0)
        a=ii;
        ip=a*Nx/wave_n; %% x(ip+1) is an integer multiple of lambda
        Np = ip+1;
        %xp = [0:Np-1]*lambda./Np;
        break
    end
end
%%
xrem = (x- lambda.*floor(x./lambda));
xtest = xrem;
[xts,idp] = sort(xtest(1:ip));
nw = Nx/ip;
lp=x(ip+1);
%%

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

convposav=zeros(Nx,Nz);
 strposav=zeros(Nx,Nz);

convnegav=zeros(Nx,Nz);
 strnegav=zeros(Nx,Nz);


mf=fullfile(baseDir,'mean_fields.mat');
 load(mf)
mp=fullfile(baseDir,'partial_sum.mat');
 load(mp)

for id = 0:nw-1%wave_n-1
	xdis  = id*lp;%id*L0;
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
	  
	  convposav = convposav + ifft((fft(convp,[],1).*kdis),[],1,'symmetric');
           strposav =  strposav + ifft((fft( strp,[],1).*kdis),[],1,'symmetric');

	  convnegav = convnegav + ifft((fft(convn,[],1).*kdis),[],1,'symmetric');
           strnegav =  strnegav + ifft((fft( strn,[],1).*kdis),[],1,'symmetric');

	  convpav = convpav + ifft((fft(convav,[],1).*kdis),[],1,'symmetric');
          strpav  =  strpav + ifft((fft( strav,[],1).*kdis),[],1,'symmetric');
end

upav = squeeze(mean(upav(1:ip,:,:),2))./nw;
vpav = squeeze(mean(vpav(1:ip,:,:),2))./nw;
wpav = squeeze(mean(wpav(1:ip,:,:),2))./nw;

oxpav =squeeze(mean(oxpav(1:ip,:,:),2))./nw;
oypav =squeeze(mean(oypav(1:ip,:,:),2))./nw;
ozpav =squeeze(mean(ozpav(1:ip,:,:),2))./nw;

viscpav =squeeze(mean(viscpav(1:ip,:,:),2))./nw;
nlpav =  squeeze(mean(  nlpav(1:ip,:,:),2))./nw;
convpav =squeeze(mean(convpav(1:ip,:,:),2))./nw;
strpav = squeeze(mean( strpav(1:ip,:,:),2))./nw;

convposav =squeeze((convposav(1:ip,:,:)))./nw;
 strposav =squeeze(( strposav(1:ip,:,:)))./nw;

convnegav =squeeze((convnegav(1:ip,:,:)))./nw;
 strnegav =squeeze(( strnegav(1:ip,:,:)))./nw; 


%%
load(fullfile(baseDir,'grid.mat'));
Xs = squeeze(X(1:ip,1,:));
Zs = squeeze(Z(1:ip,1,:));
Xrem = (Xs- lambda.*floor(Xs./lambda));

Xph  =Xrem(idp,:);
Zph  =Zs(idp,:);

uph  = upav(idp,:);
vph  = vpav(idp,:);
wph  = wpav(idp,:);
oxph = oxpav(idp,:);
oyph = oypav(idp,:);
ozph = ozpav(idp,:);
viscph=viscpav(idp,:);
nlph  =  nlpav(idp,:);
convph=convpav(idp,:);
strph = strpav(idp,:);
convposph=convposav(idp,:);
 strposph= strposav(idp,:);
convnegph=convnegav(idp,:);
 strnegph= strnegav(idp,:);


%%
close all
figure
subplot(1,2,1)
pcolor(Xs,Zs,oypav)
xlim([0 lambda])
colorbar
%shading flat

subplot(1,2,2)
pcolor(Xph,Zph,oyph)
colorbar

%%


load(fullfile(baseDir,'slines.mat'))
Zslph  =Zsl(idp,:);

load(fullfile(baseDir,'potexact.mat'))
phiph=interpolateSolution(phiexact,double(Xph),double(Zph));
phiph=reshape(phiph,[ip Nz]);

phirem = (phisl- lambda.*floor(phisl./lambda));
phislph=phirem(idp,:);

load(fullfile(baseDir,'phi_interp_2d.mat'))
uphiph=squeeze(uphi(idp,1,:));
wphiph=squeeze(wphi(idp,1,:));

%%
figure
subplot(1,2,1)
pcolor(Xph,Zph,uphiph)
%shading flat
axis equal
colorbar
ylim([-0.03 1])
clim([0.7 1.3])
subplot(1,2,2)
pcolor(Xph,Zph,wphiph)
axis equal
ylim([-0.03 1])
colorbar
clim([-0.1 0.1])

%shading flat
%shading flat
% upav = squeeze(mean(upav(1:Nxl,:,:),2))./wave_n;
% vpav = squeeze(mean(vpav(1:Nxl,:,:),2))./wave_n;
% wpav = squeeze(mean(wpav(1:Nxl,:,:),2))./wave_n;
% 
% oxpav =squeeze(mean(oxpav(1:Nxl,:,:),2))./wave_n;
% oypav =squeeze(mean(oypav(1:Nxl,:,:),2))./wave_n;
% ozpav =squeeze(mean(ozpav(1:Nxl,:,:),2))./wave_n;
% 
% viscpav =squeeze(mean(viscpav(1:Nxl,:,:),2))./wave_n;
% nlpav =  squeeze(mean(  nlpav(1:Nxl,:,:),2))./wave_n;
% convpav =  squeeze(mean(convpav(1:Nxl,:,:),2))./wave_n;
% strpav =  squeeze(mean( strpav(1:Nxl,:,:),2))./wave_n; 

 mf=fullfile(baseDir,'mean_phase_fields.mat');
 save(mf,'Xph','Zph','phiph','phislph','Zslph',...
     'uph','vph','wph','uphiph','wphiph',...
     'oxph','oyph','ozph',...
     'viscph','nlph','convph','strph',...
     'convposph','strposph','convnegph','strnegph' );
% save(mf,'upav','vpav','wpav','oxpav','oypav','ozpav','viscpav','nlpav','convpav','strpav');

%%
