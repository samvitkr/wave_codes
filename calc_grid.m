clear
close all
%baseDir = '/users/1/kuma0458/wave/wavy_ret180';
%baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';
%baseDir = '/users/1/kuma0458/wave/c2ak2_re180/run';
%baseDir = '/scratch.global/kuma0458/c24ak2_re180/run';
baseDir = '/scratch.global/kuma0458/c14ak1_re180/run';

Nx=256;
Ny=192;
Nz=128;
ak=0.1;
wave_n=12;
fn    = 'grid.h5';
fname = fullfile(baseDir,fn); 
fnmat = 'grid.mat';
fngr   = fullfile(baseDir,fnmat);
zz   = h5read(fname, '/zz');
zw   = h5read(fname, '/zw');
pex  = h5read(fname, '/pex');
pey = h5read(fname, '/pey');

k0=wave_n*pex;
a=ak/k0;
kx=pex*[0:Nx/2-1,-Nx/2:-1]';
ky=pey*[0:Ny/2-1,-Ny/2:-1]';
Lx=2*pi/pex;
Ly=2*pi/pey;
x=[0:Nx-1]*Lx/Nx;
y=[0:Ny-1]*Ly/Ny;
etaz = a*cos(k0*x)';
feta = fft(etaz);
detadx = ifft( feta.*(1i.*kx),'symmetric' );

[Zeta ,~] = meshgrid(zz,etaz);
[Zetaw ,Eta] = meshgrid(zw,etaz);
dZetadx = (Zeta-1).*( detadx./(1-etaz) );
dZetadx = reshape( dZetadx,[Nx 1 Nz] );
dZetadz = 1./(1-etaz);
zg = Eta+Zeta.*(1-Eta);
zgw = Eta+Zetaw.*(1-Eta);

Z=single(zeros(Nx,Ny,Nz));
Zw = single(zeros(Nx,Ny,Nz));
X=Z;
Y=Z;

for i =1:Nx
    for j =1:Ny
        for k=1:Nz
            X(i,j,k)=x(i);
            Y(i,j,k)=y(j);
            Z(i,j,k)=zg(i,k);
            Zw(i,j,k)=zgw(i,k);
        end
    end
end
save(fngr,'X','Y','Z','Zw','dZetadz','dZetadx','a','k0')
%%
%tstart=3020200000;
%step  =    200000;
%tend  =3052000000;

tstart=3550000000;
step = 400000;
  tend=3578000000;
for tstep=tstart:step:tend
	fn=sprintf('Sol%014d.h5',tstep)
	fname   = fullfile(baseDir,fn);   % <-- edit if naming differs
	fnmat = sprintf('grid%014d.mat',tstep);
	fng = fullfile(baseDir,fnmat);
	fprintf('Reading %s\n', fname);
	eta = h5read(fname, '/eta');
	
	%% wave numbers and fft
	etaz=eta(:,1);
	%plot(etaz)
	feta = fft(etaz);
	detadx = ifft( feta.*(1i.*kx),'symmetric' );
	
	[Zeta Eta] = meshgrid(zz,etaz);
	[Zetaw Eta] = meshgrid(zw,etaz);
	
	dZetadx = (Zeta-1).*( detadx./(1-etaz) );
	dZetadx = reshape( dZetadx,[Nx 1 Nz] );
	dZetadz = 1./(1-etaz);
	
	zg = Eta+Zeta.*(1-Eta);
	zgw = Eta+Zetaw.*(1-Eta);
	
	Z=single(zeros(Nx,Ny,Nz));
	Zw = single(zeros(Nx,Ny,Nz));
	
	for i =1:Nx
	    for j =1:Ny
	        for k=1:Nz
	            Z(i,j,k)=zg(i,k);
	            Zw(i,j,k)=zgw(i,k);
	        end
	    end
	end
	%save('grid.mat','X','Y','Z','Zw','dZetadz','dZetadx')
	save(fng,'Z','Zw','dZetadz','dZetadx')
end
