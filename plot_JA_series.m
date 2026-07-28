clear
close all
baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
%baseDir = '/scratch.global/kuma0458/c14ak1_re180/run';
load(fullfile(baseDir,'JAseries_statwave.mat'))


Nx=256;
Ny=192;
Nz=128;
fname = fullfile(baseDir,'grid.h5');
zz   = h5read(fname, '/zz');
zw   = h5read(fname, '/zw');
pex  = h5read(fname, '/pex');
pey = h5read(fname, '/pey');
Ly=2*pi/pey;
Lx=2*pi/pex;
kx=pex*[0:Nx/2-1,-Nx/2:-1]';
load(fullfile(baseDir,'grid.mat'))
load(fullfile(baseDir,'phi_interp_2d.mat'),'uphi','wphi')
	Jacobian=1./dZetadz;
	vol = Lx;

figure
subplot(1,2,1)
hold on
plot(phidots,'-r')
plot(Ts,'-b')
plot(phidots+Ts,'-k')
hold off
yline(vol)
subplot(1,2,2)
hold on
%plot(Tconvs,':b')
%plot(Tstrs,'--b')
plot(Tnls,'.-b')
plot(Tviscs,'-.r')
plot(Ts,'-.k')
hold off


%%

%trapz(X(:,1,1),Jacobian);

    uphi=squeeze(uphi);
    uphiin=trapz(zz,uphi,2);
    Jstar=trapz(X(:,1,1),Jacobian.*uphiin);
    Jstar = Jstar./vol;

    figure
plot(check./Jstar)

load(fullfile(baseDir,'flowrate.mat'))
figure
plot(flowrate)