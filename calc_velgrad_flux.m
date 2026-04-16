%% Read timestep 10 from OCF data and save into a .mat file
clear
close all
%baseDir = '/users/1/kuma0458/wave/wavy_wall';
baseDir = '/users/1/kuma0458/wave/wavy_ret180';
%baseDir = '/users/1/kuma0458/wave/wave_c_2';

%load('grid.mat')
tic
ret=180;
%ret=10;
nu=1/ret;
tstart=20200000;
step=200000;
tend=32000000;

for tstep=tstart:step:tend
    % fn=		sprintf('DAT000%03d99999999',tstep)
    % fnmat= sprintf('gradflux000%03d99999999.mat',tstep)
    %fng = sprintf('grid%014d.mat',tstep);
    %fng=fullfile(baseDir,fng);
    fng=fullfile(baseDir,'grid.mat')
    load(fng);
    fn=sprintf('Sol%014d.h5',tstep);
    fnmat=sprintf('gradflux%014d.mat',tstep);
    fname = fullfile(baseDir,fn);
    fnamemat=fullfile(baseDir,fnmat);
    % info=h5info(fname)
    fprintf('Reading %s\n', fname);

    % --- Read datasets (no leading slash) ---
    zz   = h5read(fname, '/zz');
    zw   = h5read(fname, '/zw');
    u    = h5read(fname, '/u');
    v    = h5read(fname, '/v');
    w    = h5read(fname, '/w');
    pp   = h5read(fname, '/pp');
    pex  = h5read(fname, '/pex');
    pey  = h5read(fname, '/pey');

    %time = h5read(fname, '/time');
    %dz   = h5read(fname, '/dz');
    %dzw   = h5read(fname, '/dzw');
    toc
    %% interpolate wc to cenn centers
    wc=permute(w,[3 1 2]);
    wc=interp1(zw(1:end-1),wc(1:end-1,:,:),zz);
    wc=single(permute(wc,[2 3 1]));
    %toc
    %%% wave numbers and fft
    [Nx,Ny,Nz] = size(u);
    kx=pex*[0:Nx/2-1,0,-Nx/2+1:-1]';
    ky=pey*[0:Ny/2-1,0,-Ny/2+1:-1]';

    fyu     =fft(permute(u ,[2 1 3]),[],1);
    fyv     =fft(permute(v ,[2 1 3]),[],1);
    fyw     =fft(permute(wc,[2 1 3]),[],1);
    dudy    =single(permute(ifft( fyu.*(1i.*ky),[],1,'symmetric'),[2 1 3]))     ;
    dvdy    =single(permute(ifft( fyv.*(1i.*ky),[],1,'symmetric'),[2 1 3]) )    ;
    dwdy    =single(permute(ifft( fyw.*(1i.*ky),[],1,'symmetric'),[2 1 3]) )    ;
    d2udy2  =single(permute(ifft( fyu.*(-(ky).^2),[],1,'symmetric'),[2 1 3]));
    d2vdy2  =single(permute(ifft( fyv.*(-(ky).^2),[],1,'symmetric'),[2 1 3]));
    d2wdy2  =single(permute(ifft( fyw.*(-(ky).^2),[],1,'symmetric'),[2 1 3]));

    clear fyu fyv fyw

%    dudzeta=zeros(size(u));
%    dvdzeta=zeros(size(u));
%    dwdzeta=zeros(size(u));
%    dwdzeta(:,:,2:end) = permute(diff(permute(w,[3 2 1]),1,1)./diff(zw),[3 2 1]);
%    dz=diff(zz);
%    dudzeta(:,:,1)=(u(:,:,2)-u(:,:,1))./dz(1);
%    dvdzeta(:,:,1)=(v(:,:,2)-v(:,:,1))./dz(1);
%    nz=Nz;
%    i      = 1:(nz-2);
%    dz_i   = dz(i);
%    dz_ip1 = dz(i+1);
%    sumdz  = dz_i + dz_ip1;
%    a  = 1./dz_ip1 - 1./sumdz;
%    b  = 1./dz_i   - 1./dz_ip1;
%    c  = 1./sumdz  - 1./dz_i;
%    % reshape for implicit expansion along 3rd dim
%    a  = reshape(a, 1, 1, []);
%    b  = reshape(b, 1, 1, []);
%    c  = reshape(c, 1, 1, []);
%    % Apply to u, v
%    dudzeta(:,:,2:nz-1) =a.*u(:,:,3:nz)+b.*u(:,:,2:nz-1)+c.* u(:,:,1:nz-2);
%    dvdzeta(:,:,2:nz-1) =a.*v(:,:,3:nz)+b.*v(:,:,2:nz-1)+c.* v(:,:,1:nz-2);
%    dudzeta(:,:,end)=(u(:,:,end)-u(:,:,end-1))./dz(end);
%    dvdzeta(:,:,end)=(v(:,:,end)-v(:,:,end-1))./dz(end);
%
%    dudz = single(dZetadz.*dudzeta);
%    dvdz = single(dZetadz.*dvdzeta);
%    dwdz = single(dZetadz.*dwdzeta);
%
%    d2udz2 = dudz.*0;
%    d2vdz2 = dudz.*0;
%    d2wdz2 = dudz.*0;
%
%    d2udz2(:,:,1)=(dudz(:,:,2)-dudz(:,:,1))./dz(1);
%    d2vdz2(:,:,1)=(dvdz(:,:,2)-dvdz(:,:,1))./dz(1);
%    d2wdz2(:,:,1)=(dwdz(:,:,2)-dwdz(:,:,1))./dz(1);
%
%    d2udz2(:,:,2:nz-1) =a.*dudz(:,:,3:nz)+b.*dudz(:,:,2:nz-1)+c.* dudz(:,:,1:nz-2);
%    d2vdz2(:,:,2:nz-1) =a.*dvdz(:,:,3:nz)+b.*dvdz(:,:,2:nz-1)+c.* dvdz(:,:,1:nz-2);
%    d2wdz2(:,:,2:nz-1) =a.*dwdz(:,:,3:nz)+b.*dwdz(:,:,2:nz-1)+c.* dwdz(:,:,1:nz-2);
%
%    d2udz2(:,:,end)=(dudz(:,:,end)-dudz(:,:,end-1))./dz(end);
%    d2vdz2(:,:,end)=(dvdz(:,:,end)-dvdz(:,:,end-1))./dz(end);
%    d2wdz2(:,:,end)=(dwdz(:,:,end)-dwdz(:,:,end-1))./dz(end);
%
%
%    fxu =fft(u,[],1);
%    fxv =fft(v,[],1);
%    fxw =fft(wc,[],1);
%    dudxi=ifft( fxu.*(1i.*kx),[],1,'symmetric');
%    dvdxi=ifft( fxv.*(1i.*kx),[],1,'symmetric');
%    dwdxi=ifft( fxw.*(1i.*kx),[],1,'symmetric');
%
%    dudx=single(dudxi+dZetadx.*(dudzeta));
%    dvdx=single(dvdxi+dZetadx.*(dvdzeta));
%    dwdx=single(dwdxi+dZetadx.*(dwdzeta));
%
%    fxdudx =fft(dudx,[],1);
%    fxdvdx =fft(dvdx,[],1);
%    fxdwdx =fft(dwdx,[],1);
%    d2udxi2=ifft( fxdudx.*(1i.*kx),[],1,'symmetric');
%    d2vdxi2=ifft( fxdvdx.*(1i.*kx),[],1,'symmetric');
%    d2wdxi2=ifft( fxdwdx.*(1i.*kx),[],1,'symmetric');
%
%
%    d2udx2=d2udxi2+dZetadx.*(d2udz2);
%    d2vdx2=d2vdxi2+dZetadx.*(d2vdz2);
%    d2wdx2=d2wdxi2+dZetadx.*(d2wdz2);
%    d2udz2 = dZetadz.*d2udz2;
%    d2vdz2 = dZetadz.*d2vdz2;
%    d2wdz2 = dZetadz.*d2wdz2;


% 1. FIX: Calculate vertical gradient of w, ensuring the bottom cell is included
    dwdzeta = zeros(size(u));
    dudzeta = dwdzeta;
    dvdzeta = dwdzeta;
    dwdzeta(:,:,2:end) = permute(diff(permute(w,[3 2 1]),1,1)./diff(zw),[3 2 1]);
    dwdzeta(:,:,1) = (w(:,:,2) - w(:,:,1)) ./ (zw(2) - zw(1)); % Added bottom boundary

    dz = diff(zz);
    dudzeta(:,:,1) = (u(:,:,2) - u(:,:,1)) ./ dz(1);
    dvdzeta(:,:,1) = (v(:,:,2) - v(:,:,1)) ./ dz(1);
    
    nz = Nz;
    i      = 1:(nz-2);
    dz_i   = dz(i);
    dz_ip1 = dz(i+1);
    sumdz  = dz_i + dz_ip1;
    
    a = 1./dz_ip1 - 1./sumdz;
    b = 1./dz_i   - 1./dz_ip1;
    c = 1./sumdz  - 1./dz_i;
    
    a = reshape(a, 1, 1, []);
    b = reshape(b, 1, 1, []);
    c = reshape(c, 1, 1, []);
    
    dudzeta(:,:,2:nz-1) = a.*u(:,:,3:nz) + b.*u(:,:,2:nz-1) + c.*u(:,:,1:nz-2);
    dvdzeta(:,:,2:nz-1) = a.*v(:,:,3:nz) + b.*v(:,:,2:nz-1) + c.*v(:,:,1:nz-2);
    dudzeta(:,:,end) = (u(:,:,end) - u(:,:,end-1)) ./ dz(end);
    dvdzeta(:,:,end) = (v(:,:,end) - v(:,:,end-1)) ./ dz(end);

    dudz = single(dZetadz.*dudzeta);
    dvdz = single(dZetadz.*dvdzeta);
    dwdz = single(dZetadz.*dwdzeta);

    % --- Calculate First Derivatives in x ---
    fxu = fft(u,[],1);
    fxv = fft(v,[],1);
    fxw = fft(wc,[],1);
    
    dudxi = ifft( fxu.*(1i.*kx),[],1,'symmetric');
    dvdxi = ifft( fxv.*(1i.*kx),[],1,'symmetric');
    dwdxi = ifft( fxw.*(1i.*kx),[],1,'symmetric');

    dudx = single(dudxi + dZetadx.*dudzeta);
    dvdx = single(dvdxi + dZetadx.*dvdzeta);
    dwdx = single(dwdxi + dZetadx.*dwdzeta);

    % --- Calculate Second Derivatives in x (The Chain Rule Fix) ---
    fxdudx = fft(dudx,[],1);
    fxdvdx = fft(dvdx,[],1);
    fxdwdx = fft(dwdx,[],1);
    
    d2udxi2 = ifft( fxdudx.*(1i.*kx),[],1,'symmetric');
    d2vdxi2 = ifft( fxdvdx.*(1i.*kx),[],1,'symmetric');
    d2wdxi2 = ifft( fxdwdx.*(1i.*kx),[],1,'symmetric');

    % 2. FIX: We must differentiate the first derivatives (dudx) with respect to zeta
    dudx_dzeta = zeros(size(u));
    dvdx_dzeta = zeros(size(u));
    dwdx_dzeta = zeros(size(u));

    % Bottom boundary
    dudx_dzeta(:,:,1) = (dudx(:,:,2) - dudx(:,:,1)) ./ dz(1);
    dvdx_dzeta(:,:,1) = (dvdx(:,:,2) - dvdx(:,:,1)) ./ dz(1);
    dwdx_dzeta(:,:,1) = (dwdx(:,:,2) - dwdx(:,:,1)) ./ dz(1);

    % Interior
    dudx_dzeta(:,:,2:nz-1) = a.*dudx(:,:,3:nz) + b.*dudx(:,:,2:nz-1) + c.*dudx(:,:,1:nz-2);
    dvdx_dzeta(:,:,2:nz-1) = a.*dvdx(:,:,3:nz) + b.*dvdx(:,:,2:nz-1) + c.*dvdx(:,:,1:nz-2);
    dwdx_dzeta(:,:,2:nz-1) = a.*dwdx(:,:,3:nz) + b.*dwdx(:,:,2:nz-1) + c.*dwdx(:,:,1:nz-2);

    % Top boundary
    dudx_dzeta(:,:,end) = (dudx(:,:,end) - dudx(:,:,end-1)) ./ dz(end);
    dvdx_dzeta(:,:,end) = (dvdx(:,:,end) - dvdx(:,:,end-1)) ./ dz(end);
    dwdx_dzeta(:,:,end) = (dwdx(:,:,end) - dwdx(:,:,end-1)) ./ dz(end);

    % Final analytical chain rule for the second x-derivative
    d2udx2 = d2udxi2 + dZetadx .* dudx_dzeta;
    d2vdx2 = d2vdxi2 + dZetadx .* dvdx_dzeta;
    d2wdx2 = d2wdxi2 + dZetadx .* dwdx_dzeta;

    % --- Calculate Second Derivatives in z ---
    d2udz2 = dudz.*0;
    d2vdz2 = dudz.*0;
    d2wdz2 = dudz.*0;

    d2udz2(:,:,1) = (dudz(:,:,2) - dudz(:,:,1)) ./ dz(1);
    d2vdz2(:,:,1) = (dvdz(:,:,2) - dvdz(:,:,1)) ./ dz(1);
    d2wdz2(:,:,1) = (dwdz(:,:,2) - dwdz(:,:,1)) ./ dz(1);

    d2udz2(:,:,2:nz-1) = a.*dudz(:,:,3:nz) + b.*dudz(:,:,2:nz-1) + c.*dudz(:,:,1:nz-2);
    d2vdz2(:,:,2:nz-1) = a.*dvdz(:,:,3:nz) + b.*dvdz(:,:,2:nz-1) + c.*dvdz(:,:,1:nz-2);
    d2wdz2(:,:,2:nz-1) = a.*dwdz(:,:,3:nz) + b.*dwdz(:,:,2:nz-1) + c.*dwdz(:,:,1:nz-2);

    d2udz2(:,:,end) = (dudz(:,:,end) - dudz(:,:,end-1)) ./ dz(end);
    d2vdz2(:,:,end) = (dvdz(:,:,end) - dvdz(:,:,end-1)) ./ dz(end);
    d2wdz2(:,:,end) = (dwdz(:,:,end) - dwdz(:,:,end-1)) ./ dz(end);
    
    d2udz2 = dZetadz.*d2udz2;
    d2vdz2 = dZetadz.*d2vdz2;
    d2wdz2 = dZetadz.*d2wdz2;

    viscu=single(nu.*(d2udx2+d2udy2+d2udz2));
    viscv=single(nu.*(d2vdx2+d2vdy2+d2vdz2));
    viscw=single(nu.*(d2wdx2+d2wdy2+d2wdz2));

    % % % % ox = dwdy-dvdz;
    % % % % oy = dudz-dwdx;
    % % % % oz = dvdx-dudy;
    % % % %
    % % % % voz = single(v.*oz);
    % % % % woy = single(wc.*oy);
    % % % % uoy = single(u.*oy);
    % % % % vox = single(v.*ox);
    % % %
    outFile=fnamemat;
    save(outFile,'zz','wc',...
        'dudx','dvdx','dwdx',...
        'dudy','dvdy','dwdy',...
        'dudz','dvdz','dwdz',...
        ... 'voz','woy','uoy','vox',...
        'viscu','viscv','viscw','-v7.3')

end

%%
div=dudx+dvdy+dwdz;
div2 = dudx+dwdz;
cl=20;
close all
figure
subplot(4,1,1)
pcolor(squeeze(X(:,10,:)),squeeze(Z(:,10,:)),squeeze(dudx(:,10,:) ))
shading flat
axis equal
xlim([X(1,1,1),X(end,end,end)])
ylim([min(Z,[],'all') max(Z,[],'all')])
c=colorbar;
clim([-cl cl])
ylabel('z')
ylabel(c,'du/dx')

subplot(4,1,2)
pcolor(squeeze(X(:,10,:)),squeeze(Z(:,10,:)),squeeze(squeeze(dvdy(:,10,:))))
shading flat
axis equal
xlim([X(1,1,1),X(end,end,end)])
ylim([min(Z,[],'all') max(Z,[],'all')])
c1=colorbar;
 clim([-cl cl])
ylabel('z')
ylabel(c1,'dv/dy')

subplot(4,1,3)
pcolor(squeeze(X(:,10,:)),squeeze(Z(:,10,:)),squeeze(squeeze(dwdz(:,10,:))))
shading flat
axis equal
xlim([X(1,1,1),X(end,end,end)])
ylim([min(Z,[],'all') max(Z,[],'all')])
c2=colorbar
clim([-cl cl])
ylabel('z')
ylabel(c2,'dw/dz')

subplot(4,1,4)
pcolor(squeeze(X(:,10,:)),squeeze(Z(:,10,:)),squeeze(div(:,10,:)))
shading flat
axis equal
xlim([X(1,1,1),X(end,end,end)])
ylim([min(Z,[],'all') max(Z,[],'all')])
c3=colorbar
clim([-cl cl])
ylabel('z')
ylabel(c3,'du/dx+dv/dy+dw/dz')
xlabel('x')

figure
pcolor(squeeze(X(:,:,50)),squeeze(Y(:,:,50)),squeeze(div(:,:,50)))
shading flat
axis equal
xlim([X(1,1,1),X(end,end,end)])
ylim([min(Z,[],'all') max(Z,[],'all')])
c3=colorbar
clim([-cl cl])
ylabel('z')
ylabel(c3,'du/dx+dv/dy+dw/dz')
xlabel('x')
