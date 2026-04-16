clear
baseDir = '/users/1/kuma0458/wave/wave_c_2';
% etaname=fullfile(baseDir,'eta.mat');
% load(etaname)
% etaz=eta(1:128,3);
% etax=eta(1:128,1);
% etay=eta(1:128:end,2);
% tic
% ret=10;
% nu=1/ret;

tstart=9000000;
step=200000;
tend=16000000;
figure
hold on
for tstep=tstart:step:tend
    
    pause(5)
fn=sprintf('Sol%014d.h5',tstep)
fname   = fullfile(baseDir,fn);
eta= h5read(fname, '/eta');
u= h5read(fname, '/u');
pex  = h5read(fname, '/pex');
pey = h5read(fname, '/pey');
[Nx,Ny,Nz] = size(u);
kx=pex*[0:Nx/2-1,0,-Nx/2+1:-1]';
ky=pey*[0:Ny/2-1,0,-Ny/2+1:-1]';
Lx=2*pi/pex;
Ly=2*pi/pey;

x=[0:Nx-1]*Lx/Nx;
y=[0:Ny-1]*Ly/Ny;

plot(x,eta)
pause(5)
clf
end
hold off

