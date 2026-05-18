clear
close all
baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';

tic
tstart=20200000;
step=  200000;
tend=  32000000;
Nx=128;
Ny=32;
Nz=128;
c=2;
me = [];
counter=0;

 phinl=zeros(Nx,Ny,Nz);
 phiad=zeros(Nx,Ny,Nz);
phist=zeros(Nx,Ny,Nz) ;

