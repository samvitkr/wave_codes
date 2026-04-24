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
counter=0;

 phinl=zeros(Nx,Ny,Nz);
 nla = zeros(Nx,Ny,Nz);
for tstep=tstart:tstart%step:tend
    counter=counter+1;
    fn=sprintf('Sol%014d.h5',tstep);
    fnmat=sprintf('gradflux%014d.mat',tstep);
    fname = fullfile(baseDir,fn);

    fnamemat=fullfile(baseDir,fnmat);
    fnmg = sprintf('grid%014d.mat',tstep);
    fng = fullfile(baseDir,fnmg);
    load(fng);
    Jacobian=1./dZetadz;


    fnp=sprintf('potvel%014d.mat',tstep);
    fnp=fullfile(baseDir,fnp);
    load(fnp)

    up=sqrt(uphi.^2+wphi.^2);
    uph=uphi./up;
    wph=wphi./up;

    W=sqrt(up.*Jacobian);

    fprintf('Reading %s\n', fname);
    mt=matfile(fnamemat);

    u    = h5read(fname, '/u')-c;
    v    = h5read(fname, '/v');
    w=mt.wc;
    ua = W.*(u.*uph + w.*wph);
    ub = W.*v;
    uc = W.*(-u.*wph + w.*uph);

    %load(fnamemat)
    ox = mt.dwdy-mt.dvdz;
    oy = mt.dudz-mt.dwdx;
    oz = mt.dvdx-mt.dudy;
    clear mt;
    oa = W.*(ox.*uph + oz.*wph);
    ob = W.*oy;
    oc = W.*(-ox.*wph + oz.*uph);
    clear ox oy oz u v w

    fub=fft2(ub);
    fuc=fft2(uc);

    fob=fft2(ob);
    foc=fft2(oc);

    phinl=phinl+(fub.*conj(foc)-fuc.*conj(fob))./(Nx*Ny);

    % % % % Tphinl=squeeze(mean(phinl,[1 2]));
    % % % %
    % % % %
    % % % % zz   = h5read(fname, '/zz');
          nla = (ub.*oc - uc.*ob);
    % % % %    JAnla = squeeze(mean(nla,2));
    % % % %    JAnla(isnan(JAnla))=0;
    % % % %
    % % % %    % JAnla = trapz(zz',JAnla,2);
    % % % %     % Tnla = trapz(X(:,1,1),JAnla);
    % % % % Tnla = mean(JAnla,1);
    % %  voz = single(v.*oz);
    % %   woy = single(w.*oy);
    % %   uoy = single((u).*oy);
    % %   vox = single(v.*ox);
    % %
      % %nl = nl+uphi.*(voz-woy)+wphi.*(uoy-vox);
    % % JAnl = squeeze(mean(nl,2));
    % % JAnl(isnan(JAnl))=0;
    % %
    % % % JAnl = trapz(zz',Jacobian.*JAnl,2);
    % % Tnl = trapz(X(:,1,1),Jacobian.*JAnl);
    % errnl=nl-nla;
    % me = [me;max(abs(errnl),[],'all')];
end
phinl=phinl./counter;

%%
% Tphinl=squeeze(mean(phinl,[1 2]));

dX=X(3,3,3)-X(2,2,2);
dY=Z(3,3,3)-Z(2,2,2);
Tnla =squeeze( (dX*dY).*sum(nla,[1 2]));
pex  = h5read(fname, '/pex');
pey  = h5read(fname, '/pey');
zz   = h5read(fname, '/zz');
kx=pex*[0:Nx/2-1,0,-Nx/2+1:-1]';dkx=kx(2)-kx(1);
ky=pey*[0:Ny/2-1,0,-Ny/2+1:-1]';dky=ky(2)-ky(1);

Tnlphi= squeeze(real((dkx*dky).*sum(phinl./(Nx*Ny),[1 2])));
%%
figure
hold on
plot(Tnla,zz,'-b')
plot(Tnlphi,zz,'-r')
hold off