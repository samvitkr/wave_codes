
clear
close all
baseDir = '/scratch.global/kuma0458/c0ak2_re180/run';
c=0;
%baseDir = '/scratch.global/kuma0458/c14ak1_re180/run';
% %c=14
% baseDir = '/scratch.global/kuma0458/c8ak1_re180/run';
% c=8;


tstart=4625000000;
step  =   5000000;
tend  =7820000000;

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
     
     uslice = squeeze(dy.*sum(uw,2));
     %udz=sum(uslice(:,2:end).*dzw,2);
     udz=trapz(zz,uslice,2);
     Jx=udz.*Jacobian;
     %J=mean(Jx);
     J=(dx.*sum(Jx));
    t=[t;time];
     flowrate = [flowrate;J];

    %%
% 
 fnmat = sprintf('grid%014d.mat',tstep);
 	load( fullfile(baseDir,fnmat),'dZetadz');
Jacobian_lab=1./dZetadz;
uslicel = squeeze(dy.*sum(u,2));
    udzl=trapz(zz,uslicel,2);
    Jxl=udzl.*Jacobian_lab;
    Jl=(dx.*sum(Jxl));%./Lx;
    %flowrate_lab = [flowrate_lab;Jl];


uslicel = squeeze(dy.*sum(u.*0+1,2));
    udzl=trapz(zz,uslicel,2);
    Jxl=udzl.*Jacobian_lab;
   vol=(dx.*sum(Jxl))


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

