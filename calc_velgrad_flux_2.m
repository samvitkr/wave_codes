%% Read timestep 10 from OCF data and save into a .mat file
clear
close all
%baseDir = '/users/1/kuma0458/wave/wavy_wall';
baseDir = '/users/1/kuma0458/wave/wavy_ret180';

%load('grid.mat')
tic
ret=180;
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

    dudzeta=zeros(size(u));
    dvdzeta=zeros(size(u));
    dwdzeta=zeros(size(u));
    dwdzeta(:,:,2:end) = permute(diff(permute(w,[3 2 1]),1,1)./diff(zw),[3 2 1]);
    dz=diff(zz);
    dudzeta(:,:,1)=(u(:,:,2)-u(:,:,1))./dz(1);
    dvdzeta(:,:,1)=(v(:,:,2)-v(:,:,1))./dz(1);
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
    dvdzeta(:,:,2:nz-1) =a.*v(:,:,3:nz)+b.*v(:,:,2:nz-1)+c.* v(:,:,1:nz-2);
    dudzeta(:,:,end)=(u(:,:,end)-u(:,:,end-1))./dz(end);
    dvdzeta(:,:,end)=(v(:,:,end)-v(:,:,end-1))./dz(end);

    % Reshape dZetadz for proper implicit expansion into 3D
    dZetadz_3D = reshape(dZetadz, [Nx, 1, 1]);
    dudz = single(dZetadz_3D.*dudzeta);
    dvdz = single(dZetadz_3D.*dvdzeta);

    fxu =fft(u,[],1);
    fxv =fft(v,[],1);
    fxw =fft(wc,[],1);
    dudxi=ifft( fxu.*(1i.*kx),[],1,'symmetric');
    dvdxi=ifft( fxv.*(1i.*kx),[],1,'symmetric');
    dwdxi=ifft( fxw.*(1i.*kx),[],1,'symmetric');

    dudx=single(dudxi+dZetadx.*(dudzeta));
    dvdx=single(dvdxi+dZetadx.*(dvdzeta));
    dwdx=single(dwdxi+dZetadx.*(dwdzeta));

    % --- CHANGE 1: Enforce divergence-free condition exactly ---
    dwdz = -(dudx + dvdy);

    % --- CHANGE 2: Exact Laplacian for Viscous Terms ---
    % Inverse Jacobian (Jinv = 1 - eta)
    Jinv = 1 ./ dZetadz_3D;

    % Viscous u
    F1_u = Jinv .* dudx;
    dF1_u_dxi = ifft( fft(F1_u, [], 1) .* (1i.*kx), [], 1, 'symmetric');
    F3_u = Jinv .* (dZetadx .* dudx + dZetadz_3D.^2 .* dudzeta);
    dF3_u_dzeta = zeros(size(u), 'single');
    dF3_u_dzeta(:,:,1) = (F3_u(:,:,2) - F3_u(:,:,1)) ./ dz(1);
    dF3_u_dzeta(:,:,2:nz-1) = a.*F3_u(:,:,3:nz) + b.*F3_u(:,:,2:nz-1) + c.*F3_u(:,:,1:nz-2);
    dF3_u_dzeta(:,:,end) = (F3_u(:,:,end) - F3_u(:,:,end-1)) ./ dz(end);
    lap_u = (dF1_u_dxi + dF3_u_dzeta) ./ Jinv + d2udy2;
    viscu = single(nu .* lap_u);

    % Viscous v
    F1_v = Jinv .* dvdx;
    dF1_v_dxi = ifft( fft(F1_v, [], 1) .* (1i.*kx), [], 1, 'symmetric');
    F3_v = Jinv .* (dZetadx .* dvdx + dZetadz_3D.^2 .* dvdzeta);
    dF3_v_dzeta = zeros(size(u), 'single');
    dF3_v_dzeta(:,:,1) = (F3_v(:,:,2) - F3_v(:,:,1)) ./ dz(1);
    dF3_v_dzeta(:,:,2:nz-1) = a.*F3_v(:,:,3:nz) + b.*F3_v(:,:,2:nz-1) + c.*F3_v(:,:,1:nz-2);
    dF3_v_dzeta(:,:,end) = (F3_v(:,:,end) - F3_v(:,:,end-1)) ./ dz(end);
    lap_v = (dF1_v_dxi + dF3_v_dzeta) ./ Jinv + d2vdy2;
    viscv = single(nu .* lap_v);

    % Viscous w
    F1_w = Jinv .* dwdx;
    dF1_w_dxi = ifft( fft(F1_w, [], 1) .* (1i.*kx), [], 1, 'symmetric');
    F3_w = Jinv .* (dZetadx .* dwdx + dZetadz_3D.^2 .* dwdzeta);
    dF3_w_dzeta = zeros(size(u), 'single');
    dF3_w_dzeta(:,:,1) = (F3_w(:,:,2) - F3_w(:,:,1)) ./ dz(1);
    dF3_w_dzeta(:,:,2:nz-1) = a.*F3_w(:,:,3:nz) + b.*F3_w(:,:,2:nz-1) + c.*F3_w(:,:,1:nz-2);
    dF3_w_dzeta(:,:,end) = (F3_w(:,:,end) - F3_w(:,:,end-1)) ./ dz(end);
    lap_w = (dF1_w_dxi + dF3_w_dzeta) ./ Jinv + d2wdy2;
    viscw = single(nu .* lap_w);

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
c2=colorbar;
clim([-cl cl])
ylabel('z')
ylabel(c2,'dw/dz')

subplot(4,1,4)
pcolor(squeeze(X(:,10,:)),squeeze(Z(:,10,:)),squeeze(div(:,10,:)))
shading flat
axis equal
xlim([X(1,1,1),X(end,end,end)])
ylim([min(Z,[],'all') max(Z,[],'all')])
c3=colorbar;
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
c3=colorbar;
clim([-cl cl])
ylabel('z')
ylabel(c3,'du/dx+dv/dy+dw/dz')
xlabel('x')
