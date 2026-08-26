
clear
close all
%baseDir = '/scratch.global/kuma0458/c0ak2_re180/run';
%c=0;
%baseDir = '/scratch.global/kuma0458/c14ak1_re180/run';
%c=14
baseDir = '/scratch.global/kuma0458/c8ak1_re180/run';
c=8;


%tstart=3025000000;
%step  =   5000000;
%tend  =4020000000;

%tstart=3020400000;
%step  =    400000;
%tend  =3084000000;

% tstart=3020200000;
% step  =    200000;
% tend  =3052000000;

%tstart=3550000000;
%step = 400000;
%tend=3578000000;

tstart= 4300000000;
  step=    1250000;
 tend = 4770000000;

%tstart=3825000000;
%step  =   5000000;
%tend  =4620000000;

load(fullfile(baseDir,'grid.mat'))

fn    = 'grid.h5';
fname = fullfile(baseDir,fn);
zz   = h5read(fname, '/zz');
zw   = h5read(fname, '/zw');
pex  = h5read(fname, '/pex');
pey = h5read(fname, '/pey');


t=[];
flowrate=[];
flowrate_lab=[];
dzw=diff(zw)';
Jacobian=1./dZetadz;
Nx=256;
dy = Y(4,4,4)-Y(3,3,3);
dx =X(3,1,1)-X(2,1,1);
kx=pex*[0:Nx/2-1,-Nx/2:-1]';
Lx=2*pi/pex;
Ly=2*pi/pey;
for tstep=tstart:step:tend
    fn=sprintf('Sol%014d.h5',tstep);
    fname = fullfile(baseDir,fn);
    fprintf('Reading %s\n', fname);
    % --- Read datasets (no leading slash) ---
    u    = h5read(fname, '/u');
    time = h5read(fname, '/time');
    ct=c*time;

     kd = exp((1i*ct).*kx);
     kdis = reshape(kd,[Nx,1,1]);
     uw = ifft( (fft(u,[],1).*kdis),[],1,'symmetric')-c;
     
     uslice = squeeze(dy.*trapz(uw,2));
     %udz=sum(uslice(:,2:end).*dzw,2);
     udz=trapz(zz,uslice,2);
     Jx=udz.*Jacobian;
     %J=mean(Jx);
     J=(dx.*trapz(Jx));
    t=[t;time];
     flowrate = [flowrate;J];

    %%
% 
 fnmat = sprintf('grid%014d.mat',tstep);
 	load( fullfile(baseDir,fnmat),'dZetadz');
Jacobian_lab=1./dZetadz;
uslicel = squeeze(dy.*trapz(u,2));
    udzl=trapz(zz,uslicel,2);
    Jxl=udzl.*Jacobian_lab;
    Jl=(dx.*trapz(Jxl));%./Lx;
    %flowrate_lab = [flowrate_lab;Jl];


uslicel = squeeze(dy.*trapz(u.*0+1,2));
    udzl=trapz(zz,uslicel,2);
    Jxl=udzl.*Jacobian_lab;
   vol=(dx.*trapz(Jxl))


flowrate_lab = [flowrate_lab;Jl/vol];


end
flowrate = (flowrate_lab-c)./(1-a);
%%
Jdot = flowrate.*0;
nf = length(flowrate);
Jdot(1)=(flowrate(2)-flowrate(1))/(t(2)-t(1));
for i=2:nf-1
    Jdot(i)=(flowrate(i+1)-flowrate(i-1))/(t(i+1)-t(i-1));
end
Jdot(end) = (flowrate(end)-flowrate(end-1))/(t(end)-t(end-1));
fj=fullfile(baseDir,'flowrate.mat');
 save(fj,'t','flowrate','Jdot')
%%

