clear
close all
%baseDir = '/users/1/kuma0458/wave/wavy_ret180';
%baseDir = '/users/1/kuma0458/wave/wave_ret180_c2';
%baseDir = '/users/1/kuma0458/wave/c2ak2_re180/run';
%baseDir = '/scratch.global/kuma0458/c2ak2_re180/run';
baseDir = '/scratch.global/kuma0458/c14ak1_re180/run';
% ==========================================
% 1. SETUP & GEOMETRY
% ==========================================
ak=0.1;
wave_n=12;
pex = 0.5;
Lx=2*pi/pex;
k0=wave_n*pex;
a=ak/k0;
%%
H = 1;
Nx = 256*6;   % Adjusted nodes in X to maintain resolution over 3*Lx span
Ny = 256*2;     % Nodes in Y
% Nx=256;
% Ny=128;
% Create a standard PDE Model
model = createpde();

xMin = 0;   % New Left Boundary
xMax = Lx;  % New Right Boundary

% 1A. Create a Rectangle Geometry first (to generate the mesh)
% Rectangle definition: [3, 4, x1, x2, x3, x4, y1, y2, y3, y4]'
R1 = [3, 4, xMin, xMax, xMax, xMin, 0, 0, H, H]';
g = decsg(R1, 'R1', ('R1')');
geometryFromEdges(model, g);

% 1B. Generate the Mesh (Structured-ish)
% Hmax dynamically scales with the new domain length
generateMesh(model, 'Hmax', (xMax - xMin)/(Nx-1), 'GeometricOrder', 'quadratic');

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
% eta = interp1(xin',etain',x_vals,'linear','extrap');
eta = a*cos(k0*x_vals);

% Calculate Sigma (0 at bottom, 1 at top) relative to the original flat box
sigma = y_vals / H;

% Apply the Warp: y_new = eta + sigma * (H - eta)
y_new = eta + sigma .* (H - eta);

% Update the nodes
Nodes(2, :) = y_new;

% ==========================================
% 3. RE-IMPORT WARPED GEOMETRY
% ==========================================
model_wavy = createpde();
geometryFromMesh(model_wavy, Nodes, Elements);

% Specify that we are solving Laplacian: div(grad(u)) = 0
specifyCoefficients(model_wavy, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', 0);

% ==========================================
% 4. IDENTIFY BOUNDARIES AUTOMATICALLY
% ==========================================
numEdges = model_wavy.Geometry.NumEdges;
BottomEdges = [];
TopEdges = [];
LeftEdges = [];
RightEdges = [];

for i = 1:numEdges
    % Get nodes on this edge
    [x_e, y_e] = getEdgeMidpoint(model_wavy, i);

    % Check geometric criteria against new limits
    if abs(x_e - xMin) < 1e-4
        LeftEdges = [LeftEdges, i];
    elseif abs(x_e - xMax) < 1e-4
        RightEdges = [RightEdges, i];
    elseif abs(y_e - H) < 1e-4
        TopEdges = [TopEdges, i];
    else
        BottomEdges = [BottomEdges, i];
    end
end

 % ==========================================
 % 5. APPLY BOUNDARY CONDITIONS & SOLVE
 % ==========================================
 % Dirichlet on Left (x=-Lx) and Right (x=2Lx)
 applyBoundaryCondition(model_wavy, 'dirichlet', 'Edge', LeftEdges, 'u', xMin);
 applyBoundaryCondition(model_wavy, 'dirichlet', 'Edge', RightEdges, 'u', xMax);
 
 % Neumann on Bottom and Top 
 applyBoundaryCondition(model_wavy, 'neumann', 'Edge', BottomEdges, 'g', 0);
 applyBoundaryCondition(model_wavy, 'neumann', 'Edge', TopEdges, 'g', 0);
 
 % Solve
 results = solvepde(model_wavy);
 phiexact=results;
 f=fullfile(baseDir,'potexact.mat');
 save(f,'phiexact','-v7.3');

% figure;
% pdeplot(model_wavy, 'XYData', results.NodalSolution, 'Colormap', 'jet');
% axis equal;
% xlim([xMin xMax]);
% ylim([-a 1]);
% title('Exact Nodal Velocity Potential (\phi)');
% xlabel('x');
% ylabel('z');

%% Helper function
function [xm, ym] = getEdgeMidpoint(model, edgeID)
    edgeNodes = findNodes(model.Mesh, 'region', 'Edge', edgeID);
    coords = model.Mesh.Nodes(:, edgeNodes);
    xm = mean(coords(1,:));
    ym = mean(coords(2,:));
end
