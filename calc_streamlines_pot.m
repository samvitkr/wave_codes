close all
clear 
% baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
Nx=256;
Ny=192;
Nz=128;
ak=0.2;
wave_n=12;

baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
load(fullfile(baseDir, 'grid.mat'));

fnp  = fullfile(baseDir, 'phi_interp_2d.mat'); % Saving as a 2D specific file
mp=matfile(fnp);

fn_pot  = fullfile(baseDir, 'potexact.mat');
load(fn_pot, 'phiexact');

fname = fullfile(baseDir,'grid.h5');
zz   = h5read(fname, '/zz');
zw   = h5read(fname, '/zw');
pex  = h5read(fname, '/pex');
pey = h5read(fname, '/pey');

k0=wave_n*pex;
a0=ak/k0;
kx=pex*[0:Nx/2-1,-Nx/2:-1]';
ky=pey*[0:Ny/2-1,-Ny/2:-1]';
L0=2*pi/k0;
Lx=2*pi/pex;
 X_slice    = double(squeeze(X( :,1, :))).';
 Z_slice    = double(squeeze(Z( :,1, :))).';
 uphi_slice = double(squeeze(mp.uphi(:,1, :))).';
 wphi_slice = double(squeeze(mp.wphi( :,1, :))).';

 delta = zeros(Nz,1);
zin =squeeze(Z_slice(:,1));
upin=squeeze(uphi_slice(:,1));

delta = cumtrapz(zin,upin);

%%
 
 
 istart=5;
istep=5;
iend=125;
startX = X_slice(istart:istep:iend,1);
startZ = Z_slice(istart:istep:iend,1);
startXb = X_slice(istart:istep:iend,end);
startZb = Z_slice(istart:istep:iend,end);
delta_sl = delta(istart:istep:iend);


%%
 Fu = scatteredInterpolant(X_slice(:), Z_slice(:), uphi_slice(:), 'linear', 'none');
 Fw = scatteredInterpolant(X_slice(:), Z_slice(:), wphi_slice(:), 'linear', 'none');
 ds = X_slice(3,3)-X_slice(2,2);
 max_steps = 100000;
% 
 num_lines = length(startX);
 slines=[];
 %slines=zeros(Nx,num_lines);
 %slinesb=zeros(Nx,num_lines);
% 
slz=[];
 for idx = 1:num_lines
    curr_X = startX(idx);
    curr_Z = startZ(idx);

    % Pre-allocate for speed
    path_X = NaN(max_steps, 1);
    path_Z = NaN(max_steps, 1);

    for k = 1:max_steps
        % Record current position
        path_X(k) = curr_X;
        path_Z(k) = curr_Z;

        % Find local velocity
        u = Fu(curr_X, curr_Z);
        w = Fw(curr_X, curr_Z);

        % Check boundaries / dead zones
        vel_mag = sqrt(u^2 + w^2);
        if isnan(u) || isnan(w) || vel_mag < 1e-6 ||curr_X>L0
            path_X(k+1:end) = [];
            path_Z(k+1:end) = [];
            break;
        end
        curr_X = curr_X +  ds;
        curr_Z = curr_Z + (w / u) * ds;
    end
slines=[slines,path_Z];
   % slines(:,idx)=path_Z;
 %   slines(:,idx) = 0.5*(slines(:,idx) +flip(slines(:,idx),1));


 end
% Wx = X_slice(10,:)';
% Wx =Wx./Wx(end);
% 
% 
% for idx = 1:1%num_lines
%     curr_X = startXb(idx);
%     curr_Z = startZb(idx);
% 
%     % Pre-allocate for speed
%     path_X = NaN(max_steps, 1);
%     path_Z = NaN(max_steps, 1);
% 
%     for k = 1:max_steps
%         % Record current position
%         path_X(k) = curr_X;
%         path_Z(k) = curr_Z;
% 
%         % Find local velocity
%         u = Fu(curr_X, curr_Z);
%         w = Fw(curr_X, curr_Z);
% 
%         % Check boundaries / dead zones
%         vel_mag = sqrt(u^2 + w^2);
%         if isnan(u) || isnan(w) || vel_mag < 1e-6
%             path_X(k+1:end) = [];
%             path_Z(k+1:end) = [];
%             break;
%         end
% 
%         % Step forward
%         % curr_X = curr_X + (u / vel_mag) * ds;
%         % curr_Z = curr_Z + (w / vel_mag) * ds;
%         curr_X = curr_X -  ds;
%         curr_Z = curr_Z - (w / u) * ds;
%     end
% 
%     % Draw the calculated streamline for this specific starting point
%      slinesb(:,idx)=path_Z;
% 
%     % plot(path_X, slines(:,idx), 'o-r', 'LineWidth', 1);
% 
%     % Optional: Draw a dot at the start point
%     % plot(startX(idx), startZ(idx), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 1);
% end
% slinesb=flip(slinesb,1);
% slnum = Wx.*slinesb+(1-Wx).*slines;
slnum=slines;
% % k0=2;
% % a0=0.1;
% 
 msl = mean(slnum, 1);
 asl = 0.5 * (max(slnum, [], 1) - min(slnum, [], 1));
 x_1D = X_slice(1, :)';
 sl = double(cos(k0 * x_1D) * asl + msl);
% z0=a0*cos(k0*x_1D);
 xi_vec   = double(squeeze(X(:, 1, 1)));
 psi_vec  = double(squeeze(Y(1, :, 1)));
 H_bar    = double(Z(1,1,end));
 zeta_vec = double(squeeze(Z(1,1,:) - Z(1,1,1))) ./ (H_bar - double(Z(1,1,1)));
 phi_target = Lx * (0:Nx-1)' / Nx;
 x_uniform_phi = zeros(Nx, num_lines);
 slq=sl.*0;
for idx = 1:num_lines

    z_line = sl(:, idx);
    phi_current = interpolateSolution(phiexact, x_1D, z_line);
    phi_current = Lx * (phi_current - phi_current(1)) / (phi_current(end) - phi_current(1));
    [phi_clean, unique_idx] = unique(phi_current);
    x_clean = x_1D(unique_idx);
    x_uniform_phi(:, idx) = interp1(phi_clean, x_clean, phi_target, 'spline', 'extrap');
    slq(:,idx)=cos(k0 * x_uniform_phi(:,idx)) * asl(idx) + msl(idx);
end
xq=x_uniform_phi;
% %xq = double(X_slice(istart:istep:iend, :)');
% %c=2;
% %%

%%


mslines=fullfile(baseDir,'slines.mat');
%save(mslines,'sl','slq','asl','msl','k0','x_1D','a0','xi_vec','psi_vec','zeta_vec','xq','x_uniform_phi','delta','delta_sl');
save(mslines,'k0','x_1D','a0','xi_vec','psi_vec','zeta_vec','xq','x_uniform_phi','delta','delta_sl');...'-append')


% %%
% hold on
% %plot(X_slice(istart:istep:iend,:)',slnum,':b')
% %plot(X_slice(istart:istep:iend,:)',sl,'-sk')
% %plot(x_1D,z0,'-k','LineWidth',1.5)
% plot(x_uniform_phi,slq,'-ok')
% hold off

