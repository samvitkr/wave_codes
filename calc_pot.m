close all
clear
m=matfile('../wave_c_2/grid.mat');
etain=m.Z(:,1,1);
[a b] = max(etain);
zstart=min(etain);
etain = circshift(etain,-b);
xin=m.X(:,1,1)-m.X(1,1,1);

xv = squeeze(m.X(1:end,1,1:end));
zv =  squeeze(m.Z(1:end,1,1:end));
% ==========================================
% 1. SETUP & GEOMETRY
% ==========================================
Lx = 2*pi;
H = 1;
Nx = 256; % Nodes in X (Matches your file resolution)
Ny = 256; % Nodes in Y
Lx=2*Lx;
Nx=2*Nx;
% Create a standard PDE Model
model = createpde();

% 1A. Create a Rectangle Geometry first (to generate the mesh)
% We start with a flat box [0, Lx] x [0, H]
R1 = [3, 4, 0, Lx, Lx, 0, 0, 0, H, H]'; 
g = decsg(R1, 'R1', ('R1')');
geometryFromEdges(model, g);

% 1B. Generate the Mesh (Structured-ish)
% We force a specific element size to get close to your grid resolution
generateMesh(model, 'Hmax', Lx/(Nx-1), 'GeometricOrder', 'quadratic');

% Extract the nodes and elements
Nodes = model.Mesh.Nodes; % 2 x NumNodes
Elements = model.Mesh.Elements; % 3 x NumTriangles

% ==========================================
% 2. WARP THE MESH (The "Transform")
% ==========================================
% Nodes(1,:) is X, Nodes(2,:) is Y
x_vals = Nodes(1, :);
y_vals = Nodes(2, :);

% Define your Bottom Curve eta(x)
% REPLACE THIS with your file reading logic: eta = interp1(file_x, file_y, x_vals);
%eta = 0.1 * sin(2 * x_vals); 

%eta = interp1(xin',etain',x_vals,'linear','extrap');
eta = 0.1*cos(2*x_vals);

% Calculate Sigma (0 at bottom, 1 at top) relative to the original flat box
% Since the mesh was created on [0,H], y_vals goes 0->H.
sigma = y_vals / H;

% Apply the Warp: y_new = eta + sigma * (H - eta)
y_new = eta + sigma .* (H - eta);

% Update the nodes
Nodes(2, :) = y_new;

% ==========================================
% 3. RE-IMPORT WARPED GEOMETRY
% ==========================================
% Crucial Step: We create a new geometry from the moved nodes.
% This allows MATLAB to recognize the new wavy boundaries.
model_wavy = createpde();
geometryFromMesh(model_wavy, Nodes, Elements);
% [FIX] Specify that we are solving Laplacian: div(grad(u)) = 0
specifyCoefficients(model_wavy, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 0);

% ==========================================
% 4. IDENTIFY BOUNDARIES AUTOMATICALLY
% ==========================================
% Since edges get re-numbered, we find them by their coordinates.
% We check the midpoint of each edge.

numEdges = model_wavy.Geometry.NumEdges;
BottomEdges = [];
TopEdges = [];
LeftEdges = [];
RightEdges = [];

for i = 1:numEdges
    % Get nodes on this edge
    [x_e, y_e] = getEdgeMidpoint(model_wavy, i);
    
    % Check geometric criteria
    % 1. Check if it's on the left/right walls
    if abs(x_e - 0) < 1e-4
        LeftEdges = [LeftEdges, i];
    elseif abs(x_e - Lx) < 1e-4
        RightEdges = [RightEdges, i];
    % 2. Check if it's Top (y ~ H)
    elseif abs(y_e - H) < 1e-4
        TopEdges = [TopEdges, i];
    % 3. Otherwise, it must be the wavy bottom
    else
        % Optional: verify it lies on eta(x)
        BottomEdges = [BottomEdges, i];
    end
end

% ==========================================
% 5. APPLY BOUNDARY CONDITIONS & SOLVE
% ==========================================

% Dirichlet on Left (u=10) and Right (u=0)
applyBoundaryCondition(model_wavy, 'dirichlet', 'Edge', LeftEdges, 'u', -1);
applyBoundaryCondition(model_wavy, 'dirichlet', 'Edge', RightEdges, 'u', 1);

% Neumann on Bottom (Flux g=5) and Top (Flux g=0)
% Note: 'g' in MATLAB is the flux term.
applyBoundaryCondition(model_wavy, 'neumann', 'Edge', BottomEdges, 'g', 0);
applyBoundaryCondition(model_wavy, 'neumann', 'Edge', TopEdges, 'g', 0);

% Solve
results = solvepde(model_wavy);
phiexact=results;
% ==========================================
% 6. VISUALIZE & EXPORT
% ==========================================
%%



% ==========================================
% 7. COMPUTE & PLOT SECOND DERIVATIVES (FIXED)
% ==========================================
% 1. Create the structured grid
x_vec = linspace(0, Lx, Nx);
y_vec = linspace(double(zstart), 1.1*H, 2*Ny);

% [CRITICAL] You must use meshgrid to create 2D matrices
[Xq, Yq] = meshgrid(x_vec, y_vec); 

% Check: Xq and Yq should be 129x128 matrices, NOT vectors.
% If they are vectors, gradient will fail.

% 2. Interpolate
u_grid = interpolateSolution(results, Xq, Yq);
Phi = interpolateSolution(results, double(xv)+pi, double(zv));
Phi = reshape(Phi,size(xv));
% [DEBUG] Force reshape just in case specific MATLAB versions flatten the outp
if isvector(u_grid)
    u_grid = reshape(u_grid, size(Xq));
end

% 3. Compute Derivatives
dx = x_vec(2) - x_vec(1);
dy = y_vec(2) - y_vec(1);

% Now gradient sees a 2D matrix and will happily return 2 outputs
[u_x, u_y] = gradient(u_grid, dx, dy);

% Second derivatives
[u_xx, u_xy] = gradient(u_x, dx, dy);
[u_yx, u_yy] = gradient(u_y, dx, dy);

%%
XG=squeeze(m.X(1:end,2,1:end))';
ZG=squeeze(m.Zw(1:end,2,1:end))';

% Interpolate velocity components from (Xq,Yq) -> (XG,ZG)
ux_grid = interp2(Xq, Yq, u_x, XG, ZG, 'linear', 0);
uy_grid = interp2(Xq, Yq, u_y, XG, ZG, 'linear', 0);
uin=ux_grid(:,1);
uin(isnan(uin))=0;
 J=trapz(ZG(:,1),uin);
uphi=ux_grid'./J;
wphi=uy_grid'./J;

x = squeeze(m.X(1:end,2,1:end));
z = squeeze(m.Zw(1:end,2,1:end));

%%
% cl = 0.1;
% close all
% figure
% subplot(4,1,1)
% pdeplot(model_wavy, 'XYData', results.NodalSolution, 'Colormap', 'jet');
% title('Solution on Wavy Domain');
% axis equal;
% xlim([0 Lx])
% ylim([min(etain) 1])
% cp=colorbar;
% ylabel(cp,'\phi')
% 
% subplot(4,1,2)
% pcolor(Xq,Yq,u_xx)
% shading flat
% axis equal
% xlim([0 Lx])
% ylim([min(etain) 1])
% clim([-cl cl])
% c=colorbar;
% ylabel(c,'du_{\phi}/dx')
% 
% subplot(4,1,3)
% pcolor(Xq,Yq,u_yy)
% shading flat
% axis equal
% xlim([0 Lx])
% ylim([min(etain) 1])
% c1=colorbar;
% clim([-cl cl])
% ylabel(c1,'dw_{\phi}/dz')
% 
% subplot(4,1,4)
% pcolor(Xq,Yq,u_xx+u_yy)
% axis equal
% shading flat
% xlim([0 Lx])
% ylim([min(etain) 1])
% c2=colorbar;
% clim([-cl cl])
% ylabel(c2,'du_{\phi}/dx + dw_{\phi}/dz')

%%
figure
subplot(2,1,1)
pcolor(x,z,uphi);
title('Solution on Wavy Domain');
axis equal;
xlim([0 Lx])
ylim([min(etain) 1])
cp=colorbar;
shading flat
ylabel(cp,'u_{\phi}')
%clim([0 0.5])

subplot(2,1,2)
pcolor(x,z,wphi);
title('Solution on Wavy Domain');
axis equal;
xlim([0 Lx])
ylim([min(etain) 1])
cp=colorbar;
shading flat
ylabel(cp,'v_{\phi}')
%clim([-0.1 0.1])

outFile='pot.mat';
    save(outFile,'Phi','uphi','wphi','-v7.3')

%%
% Helper function to find edge midpoints
function [xm, ym] = getEdgeMidpoint(model, edgeID)
    % Get nodes for this edge
    edgeNodes = findNodes(model.Mesh, 'region', 'Edge', edgeID);
    % Average their coordinates to find center
    coords = model.Mesh.Nodes(:, edgeNodes);
    xm = mean(coords(1,:));
    ym = mean(coords(2,:));
end
