
clear
close all
baseDir = '/scratch.global/kuma0458/c-2ak2_re180/run';
c=-2;
%baseDir = '/scratch.global/kuma0458/c14ak1_re180/run';
%c=14
%baseDir = '/scratch.global/kuma0458/c8ak1_re180/run';
%c=8;


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

%tstart=4021250000;
%step =    1250000;
%tend = 4270000000;

tstart=3825000000;
step =    5000000;
tend =4620000000;

load(fullfile(baseDir,'grid.mat'))

fn    = 'grid.h5';
fname = fullfile(baseDir,fn);
zz   = h5read(fname, '/zz');
zw   = h5read(fname, '/zw');
pex  = h5read(fname, '/pex');
pey = h5read(fname, '/pey');


dzw =diff(zw)';
Jacobian=1./dZetadz;
t=[];
flowrate=[];
dzw=diff(zw)';
Jacobian=1./dZetadz;
Nx=256;
dy = Y(4,4,4)-Y(3,3,3);
dx =X(3,1,1)-X(2,1,1);
pex=0.5;
kx=pex*[0:Nx/2-1,-Nx/2:-1]';

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
    u = ifft( (fft(u,[],1).*kdis),[],1,'symmetric')-c;

    uslice = squeeze(dy.*sum(u,2));
    udz=sum(uslice(:,2:end).*dzw,2);
    Jx=udz.*Jacobian;
    J=mean(Jx);
    t=[t;time];
    flowrate = [flowrate;J];
end

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

%plot(t,flowrate,'o-')
