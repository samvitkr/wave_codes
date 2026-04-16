clear
close all
load('phi.mat')
load('../wave_c_2/grid.mat')
baseDir = '/users/1/kuma0458/wave/wave_c_2';

%load('grid.mat')
tic
ret=10;
nu=1/ret;
t=9000000;
fn=sprintf('Sol%014d.h5',t);
fname = fullfile(baseDir,fn);
zz   = h5read(fname, '/zz');
pex  = h5read(fname, '/pex');
pey  = h5read(fname, '/pey');
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
uin=squeeze(uphi(1,1,:));
uin(isnan(uin))=0;
Js=trapz(squeeze(Z(1,1,:)),uin)
uphi=uphi./Js;

save('phi.mat','Phi','uphi','wphi','-v7.3')
%%
%set(h_slice,'Color','b','LineWidth',1.2)
