clear
close all
load('potexact.mat')
%load('../wave_c_2/grid.mat')
%baseDir = '/users/1/kuma0458/wave/wave_c_2_ret180';
baseDir = '/users/1/kuma0458/wave/wavy_ret180';

%cwave=2;
%load('grid.mat')
cwave=0;
tic
ret=10;
nu=1/ret;
%tstep=9000000;
tstart=20200000;
step=200000;
tend=20200000;
for tstep=tstart:step:tend
%fng = sprintf('grid%014d.mat',tstep);
fng=sprintf('grid.mat');
fng=fullfile(baseDir,fng);
   load(fng);
fn=sprintf('Sol%014d.h5',tstep)
fname = fullfile(baseDir,fn);
%fnp=sprintf('potvel%014d.mat',tstep);
fnp='potvel.mat';
fnp=fullfile(baseDir,fnp);
zz   = h5read(fname, '/zz');
    zw = h5read(fname,'/zw');
pex  = h5read(fname, '/pex');
pey  = h5read(fname, '/pey');
t = h5read(fname, '/time');

xshift=cwave*t - 2*pi*floor(cwave*t/(2*pi));
xv = squeeze(X(1:end,1,1:end));
zv =  squeeze(Z(1:end,1,1:end));
Phi = interpolateSolution(phiexact, double(xv)-xshift, double(zv));
Phi = reshape(Phi,size(xv))+xshift;

u  = Phi;
[Nx,Nz] = size(u);
u=reshape(Phi,Nx, 1, Nz)-X(:,1,:);
kx=pex*[0:Nx/2-1,0,-Nx/2+1:-1]';
dudzeta=zeros(size(u));
dz=diff(zz);
dudzeta(:,:,1)=(u(:,:,2)-u(:,:,1))./dz(1);
nz=Nz;
i      = 1:(nz-2);
dz_i   = dz(i);
dz_ip1 = dz(i+1);
sumdz  = dz_i + dz_ip1;
a  = 1./dz_ip1 - 1./sumdz;
b  = 1./dz_i   - 1./dz_ip1;
c  = 1./sumdz  - 1./dz_i;
% reshape for implicit expansion along 3rd dim
a  = reshape(a, 1, 1, []);
b  = reshape(b, 1, 1, []);
c  = reshape(c, 1, 1, []);
% Apply to u, v
dudzeta(:,:,2:nz-1) =a.*u(:,:,3:nz)+b.*u(:,:,2:nz-1)+c.* u(:,:,1:nz-2);
dudz = single(dZetadz.*dudzeta);
fxu =fft(u,[],1);
dudxi=ifft( fxu.*(1i.*kx),[],1,'symmetric');

dudx=single(dudxi+dZetadx.*(dudzeta));
uphi=dudx+1.0;
wphi=dudz;
dzw =diff(zw)';
    Jacobian=1./dZetadz;

    %     uin=squeeze(u(1,:,:));
    %     yin=squeeze(Y(1,:,:));
    %     zin=squeeze(Z(1,1,:));
         % dy = Y(4,4,4)-Y(3,3,3);
         % bot = 0.1*cos(2*(X(:,1,1)-c*time));
    %
    uslice = reshape(uphi,128,128);
    uslice(isnan(uslice))=0;
    udz=sum(uslice(:,2:end).*dzw,2);
    Jx=udz.*Jacobian;
%    xshift = c*time;
    dx =X(3,1,1)-X(2,1,1);

    %ishift=round(xshift/dx);
    J=mean(Jx);
% uin=squeeze(uphi(1,1,:));
% uin(isnan(uin))=0;
uphi=uphi./J;
Phi=Phi./J;
wphi=wphi./J;

%
 save(fnp,'Phi','uphi','wphi','-v7.3')

%%%
%set(h_slice,'Color','b','LineWidth',1.2)
end