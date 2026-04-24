
clear
%close all
%baseDir = fullfile('/projects', 'standard', 'shenl','naras062','JA-Samvit','Retau-10','RUN-c0');
%   datadir=fullfile(getenv('MSIPROJECT'),'shared','kuma0458','open_channel_flow_180','data');
%baseDir = '/users/1/kuma0458/wave/wavy_wall';
baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';

%baseDir = '/users/1/kuma0458/wave/wavy_ret180';


tstart=20200000;
step=200000;
tend=32000000;
c=2;
t=[];
flowrate=[];
for tstep=tstart:step:tend
        fng = sprintf('grid%014d.mat',tstep);
    %fng='grid.mat';
       fng = fullfile(baseDir,fng)
    load(fng);
    fn=sprintf('Sol%014d.h5',tstep);
    %fnmat=sprintf('gradflux%014d.mat',tstep);
    fname = fullfile(baseDir,fn);
    %fnamemat=fullfile(baseDir,fnmat);
    %mg=matfile(fnamemat);

    % dudx=mg.dudx;
    % dvdy=mg.dvdy;
    % dwdz=mg.dwdz;
    % div=dudx+dvdy+dwdz;
    % info=h5info(fname)
    fprintf('Reading %s\n', fname);
    % --- Read datasets (no leading slash) ---
    u    = h5read(fname, '/u');
    zz = h5read(fname,'/zz');
    zw = h5read(fname,'/zw');

    time = h5read(fname, '/time');
    dzw =diff(zw)';
    Jacobian=1./dZetadz;

    %     uin=squeeze(u(1,:,:));
    %     yin=squeeze(Y(1,:,:));
    %     zin=squeeze(Z(1,1,:));
         dy = Y(4,4,4)-Y(3,3,3);
         bot = 0.1*cos(2*(X(:,1,1)-c*time));
    %
    uslice = squeeze(dy.*sum(u,2));
    udz=sum(uslice(:,2:end).*dzw,2);
    Jx=udz.*Jacobian;
    xshift = c*time;
    dx =X(3,1,1)-X(2,1,1);

    %ishift=round(xshift/dx);
    J=mean(Jx);
     t=[t;time];
     flowrate = [flowrate;J];
end

%%
% uslice = squeeze(squeeze(sum(u,2)))*dy;
%
% jac = (1-bot);
%
% for ii =1:128
%     upj=(uslice(ii,2:end))';
%     Jx(ii)=sum(upj.*dzw);
% end
% Jx=Jx';
%Jx=jac.*Jx;
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

plot(t,flowrate,'o-')