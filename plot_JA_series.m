clear
% close all
baseDir = '/scratch.global/kuma0458/c8ak1_re180/run';
fign =                     "JAcheck_c8ak1_re180_static.fig";

%baseDir = '/scratch.global/kuma0458/c8ak1_re180/run';
load(fullfile(baseDir,'JAseries_statwave.mat'))
Nx=256;
Ny=192;
Nz=128;

% tstart=3550000000;
% step  =    400000;
% tend  =3578000000;

% tstart=3025000000;
% step  =   5000000;
% % tend  =4020000000;
% 
tstart=4300000000;
step =    1250000;
tend = 4770000000;
% % 
% tstart=3825000000;
% step =    5000000;
% tend =4620000000;

t = [tstart:step:tend].*(1e-8);

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

f=figure
subplot(1,2,1)
hold on
plot(t,phidots,'-r')
plot(t,Ts,'-.b')
plot(t,phidots+Ts,'-k')
hold off
yline(vol)
legend('(L_x/J_*)dJ/dt','T','T + (L_x/J_*)dJ/dt')
xlabel('time')
legend boxoff
subplot(1,2,2)
hold on

plot(t,check)
hold off
xlabel('time')
ylabel('%error')

saveas(f,fign);


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
%
plot(flowrate./(2*pi),'o-')