close all
clear 
Nx=128;
Nz=128;
load('grid.mat')
load('potexact.mat')
%load('../wave_c_2/grid.mat')
baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';
xv = squeeze(X(1:end,1,1:end));
zv =  squeeze(Z(1:end,1,1:end));
Phi = interpolateSolution(phiexact, double(xv), double(zv));
Phi = reshape(Phi,size(xv));
Phi = reshape(Phi,Nx, 1, Nz);
c=2;
tstep=32000000;
time = tstep/(1e+8);
ct=c*time;
pex=1;
kx=pex*[0:Nx/2-1,0,-Nx/2+1:-1]';
kd = exp((1i*ct).*kx);
kdis = reshape(kd,[Nx,1,1]);
fnp=sprintf('potvel%014d.mat',tstep);
%fnp=sprintf('potvel.mat');
    fnp=fullfile(baseDir,fnp);
    mp=matfile(fnp);
uphi = ifft( (fft(mp.uphi,[],1).*kdis),[],1,'symmetric');
wphi = ifft( (fft(mp.wphi,[],1).*kdis),[],1,'symmetric');
save('potvel.mat','Phi','uphi','wphi');
X_slice    = double(squeeze(X( :,1, :))).';
Z_slice    = double(squeeze(Z( :,1, :))).';
uphi_slice = double(squeeze(uphi(:,1, :))).';
wphi_slice = double(squeeze(wphi( :,1, :))).';
istart=9;
istep=3;
iend=125;
startX = X_slice(istart:istep:iend,1);
startZ = Z_slice(istart:istep:iend,1);
startXb = X_slice(istart:istep:iend,end);
startZb = Z_slice(istart:istep:iend,end);
Fu = scatteredInterpolant(X_slice(:), Z_slice(:), uphi_slice(:), 'linear', 'none');
Fw = scatteredInterpolant(X_slice(:), Z_slice(:), wphi_slice(:), 'linear', 'none');
domain_width = max(X_slice(:)) - min(X_slice(:));
%ds = domain_width / 500;
ds = X_slice(3,3)-X_slice(2,2);
max_steps = 10000;

num_lines = length(startX);
slines=zeros(Nx,num_lines);
slinesb=zeros(Nx,num_lines);

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
        if isnan(u) || isnan(w) || vel_mag < 1e-6
            path_X(k+1:end) = [];
            path_Z(k+1:end) = [];
            break;
        end
        curr_X = curr_X +  ds;
        curr_Z = curr_Z + (w / u) * ds;
    end

    slines(:,idx)=path_Z;
    slines(:,idx) = 0.5*(slines(:,idx) +flip(slines(:,idx),1));


end
Wx = X_slice(10,:)';
Wx =Wx./Wx(end);


for idx = 1:num_lines
    curr_X = startXb(idx);
    curr_Z = startZb(idx);

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
        if isnan(u) || isnan(w) || vel_mag < 1e-6
            path_X(k+1:end) = [];
            path_Z(k+1:end) = [];
            break;
        end

        % Step forward
        % curr_X = curr_X + (u / vel_mag) * ds;
        % curr_Z = curr_Z + (w / vel_mag) * ds;
        curr_X = curr_X -  ds;
        curr_Z = curr_Z - (w / u) * ds;
    end

    % Draw the calculated streamline for this specific starting point
    slinesb(:,idx)=path_Z;

    % plot(path_X, slines(:,idx), 'o-r', 'LineWidth', 1);

    % Optional: Draw a dot at the start point
    % plot(startX(idx), startZ(idx), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 1);
end
slinesb=flip(slinesb,1);
slnum = Wx.*slinesb+(1-Wx).*slines;
k0=2;
a0=0.1;

msl = mean(slnum, 1);
asl = 0.5 * (max(slnum, [], 1) - min(slnum, [], 1));
x_1D = X_slice(1, :)';
sl = cos(k0 * x_1D) * asl + msl;
z0=a0*cos(k0*x_1D);
xi_vec   = double(squeeze(X(:, 1, 1)));
psi_vec  = double(squeeze(Y(1, :, 1)));
H_bar    = double(Z(1,1,end));
zeta_vec = double(squeeze(Z(1,1,:) - Z(1,1,1))) ./ (H_bar - double(Z(1,1,1)));
xq = double(X_slice(istart:istep:iend, :)');
c=2;
%%
save('slines.mat','sl','asl','msl','k0','x_1D','a0','xi_vec','psi_vec','zeta_vec','xq','c');

hold on
plot(X_slice(istart:istep:iend,:)',slnum,':b')
plot(X_slice(istart:istep:iend,:)',sl,'-k')
plot(x_1D,z0,'-k','LineWidth',1.5)
hold off

