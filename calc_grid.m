%% Read timestep 10 from OCF data and save into a .mat file
clear
close all
baseDir = '/users/1/kuma0458/wave/wavy_ret180';

tstart=20200000;
step  =20200000;
tend  =20200000;
figure
hold on
for tstep=tstart:tstart...step:tend
% fn=sprintf('DAT000%03d99999999',tstep)
fn=sprintf('Sol%014d.h5',tstep)
fname   = fullfile(baseDir,fn);   % <-- edit if naming differs
% info=h5info(fname)
fnmat = sprintf('grid%014d.mat',tstep);
%fng = fullfile(baseDir,fnmat);
fng=fullfile(baseDir,'grid.mat')
fprintf('Reading %s\n', fname);

% %     % --- Read datasets (no leading slash) ---
zz   = h5read(fname, '/zz');
zw   = h5read(fname, '/zw');
u    = h5read(fname, '/u');
w    = h5read(fname, '/w');

pex  = h5read(fname, '/pex');
pey = h5read(fname, '/pey');
eta = h5read(fname, '/eta');
[Nx,Ny,Nz] = size(u);


%% wave numbers and fft
[Nx,Ny,Nz] = size(u);
kx=pex*[0:Nx/2-1,0,-Nx/2+1:-1]';
ky=pey*[0:Ny/2-1,0,-Ny/2+1:-1]';
Lx=2*pi/pex;
Ly=2*pi/pey;

x=[0:Nx-1]*Lx/Nx;
y=[0:Ny-1]*Ly/Ny;
etaz=eta(:,1);
plot(etaz)
feta = fft(etaz);
detadx = ifft( feta.*(1i.*kx),'symmetric' );

[Zeta Eta] = meshgrid(zz,etaz);
[Zetaw Eta] = meshgrid(zw,etaz);

dZetadx = (Zeta-1).*( detadx./(1-etaz) );
dZetadx = reshape( dZetadx,[Nx 1 Nz] );
dZetadz = 1./(1-etaz);

zg = Eta+Zeta.*(1-Eta);
zgw = Eta+Zetaw.*(1-Eta);
X=single(zeros(Nx,Ny,Nz));
Y=single(zeros(Nx,Ny,Nz));
Z=single(zeros(Nx,Ny,Nz));
Zw = single(zeros(Nx,Ny,Nz));

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
%save('grid.mat','X','Y','Z','Zw','dZetadz','dZetadx')
 save(fng,'X','Y','Z','Zw','dZetadz','dZetadx')
end
