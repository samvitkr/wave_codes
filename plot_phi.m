
baseDir = '/scratch.global/kuma0458/c8ak1_re180/run';
load(fullfile(baseDir,'grid.mat'))
load(fullfile(baseDir,'phi_interp_2d.mat'))
load(fullfile(baseDir,'slines.mat'))

close all

Xs=squeeze(X(:,1,:));
Zs=squeeze(Z(:,1,:));



fn=figure;
hold on
contourf(Xs,Zs,Phi_2D,[0:pi/60:pi],'LineWidth',1,'Color','w')
plot(x_uniform_phi,slq,'-k','LineWidth',1)
hold off
shading interp
colormap turbo
axis equal

% --- Y-axis limits from original script ---

% --- New X-axis formatting ---
% Set limits from 0 to pi
xlim([0, pi])

% Set tick marks at multiples of pi/4
xticks([0, pi/6, pi/3, pi/2, 2*pi/3, 5*pi/6, pi])
%xticks(0:pi/12:pi)

% Set the labels to the unsimplified pi/12 fractions

% Set the labels using TeX formatting
xticklabels({'0', '\pi/6', '2\pi/6', '3\pi/6', '4\pi/6', '5\pi/6', '\pi'})

% Ensure the labels render as TeX 
set(gca, 'TickLabelInterpreter', 'tex')

% Set tick direction to point outward
set(gca, 'TickDir', 'out')
clim([0 pi])
set(gca,'FontSize',12)
ylim([-a0 1])

cb=colorbar;
ylabel(cb,'$\varphi/H$','Interpreter','latex','FontSize',14)
cb.Ticks = 0:pi/6:pi; % Set the tick locations at intervals of pi/6

% Set the simplified pi/6 fraction labels
cb.TickLabels = {'0', '\pi/6', '2\pi/6', '3\pi/6', '4\pi/6', '5\pi/6', '\pi'};

% Ensure the colorbar labels also render as TeX
cb.TickLabelInterpreter = 'tex';
xlabel('x/H')
ylabel('z/H')
% % close all
% % clear
% % L0=pi/3;
% % xt=[0:0.25:1].*L0;
% % a0=0.033;


saveas(fn,fullfile(baseDir,'phi.fig'))

