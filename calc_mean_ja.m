clear
close all
baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';

Nx=128;
Ny=32;
Nz=128;
nlav=zeros(Nx,Ny,Nz);
viscav=zeros(Nx,Ny,Nz);
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
    kdis = exp((1i*ct).*kx);
	% fn=               sprintf('DAT000%03d99999999',tstep)
    % fnmat= sprintf('gradflux000%03d99999999.mat',tstep)
    fngs = sprintf('jafields%014d.mat',tstep)
    fng = fullfile(baseDir,fngs);
    load(fng)
	fvisc = fft(JAvisc,[],1).*kdis;
	viscav = viscav+ifft(fvisc,[],1,'symmetric');
	clear fvisc
	fnl = fft(JAnl,[],1).*kdis;
        nlav = nlav+ifft(fnl,[],1,'symmetric');
       	clear fnl
	counter =counter+1;
end
nlav=nlav./counter;
viscav=viscav./counter;
save('mean_JA.mat','nlav','viscav');
