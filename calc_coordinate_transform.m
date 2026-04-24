clear
close all
% baseDir = fullfile('/projects', 'standard', 'shenl','naras062','JA-Samvit','Retau-10','RUN-c0');
%baseDir = '/users/1/kuma0458/wave/wavy_wall';
%baseDir = '/users/1/kuma0458/wave/wave_c_2';
%baseDir = '/users/1/kuma0458/wave/wavy_ret180';
baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';

%load('grid.mat')
tic
%ret=10;
ret=180;
nu=1/ret;
tstart=20200000;
step=  200000;
tend=  32000000;
% ff=fullfile(baseDir,'phi.mat');
fj=fullfile(baseDir,'flowrate.mat')
load(fj);
Ts = Jdot.*0;
phidots=Jdot.*0;
Tnls=Ts;
Tviscs=Ts;
check=Jdot.*0;
Nx=128;
Ny=32;
Nz=128;
%uphi=reshape(uphi,Nx,1,Nz);
%wphi=reshape(wphi,Nx,1,Nz);
c=2;
me = [];
phinl=zeros(Nx,Ny,Nz);
for tstep=tstart:step:tend
    fn=sprintf('Sol%014d.h5',tstep);
    fnmat=sprintf('gradflux%014d.mat',tstep);
    fname = fullfile(baseDir,fn);

    fnamemat=fullfile(baseDir,fnmat);
    fnmg = sprintf('grid%014d.mat',tstep);
    fng = fullfile(baseDir,fnmg);
    load(fng);
    fnp=sprintf('potvel%014d.mat',tstep);
    fnp=fullfile(baseDir,fnp);
    load(fnp)

    up=sqrt(uphi.^2+wphi.^2);
    uph=uphi./up;
    wph=wphi./up;

    fprintf('Reading %s\n', fname);
    mt=matfile(fnamemat);

    u    = h5read(fname, '/u')-c;
    v    = h5read(fname, '/v');
    w=mt.wc;
    ua = up.*(u.*uph + w.*wph);
    ub = up.*v;
    uc =up.*(-u.*wph + w.*uph);
    
    %load(fnamemat)
    ox = mt.dwdy-mt.dvdz;
    oy = mt.dudz-mt.dwdx;
    oz = mt.dvdx-mt.dudy;
    clear mt;
    oa = ox.*uph + oz.*wph;
    ob = oy;
    oc =-ox.*wph + oz.*uph;

    fub=fft2(ub);
    fuc=fft2(uc);

    fob=fft2(ob);
    foc=fft2(oc);

    phinl=phinl+(fub.*conj(foc)-fuc.*conj(fub))./(Nx*Ny);


    % nla = (ub.*oc - uc.*ob);
%     voz = single(v.*oz);
%     woy = single(w.*oy);
%     uoy = single((u).*oy);
%     vox = single(v.*ox);
%     nl = uphi.*(voz-woy)+wphi.*(uoy-vox);
% errnl=nl-nla;
% me = [me;max(abs(errnl),[],'all')];
end