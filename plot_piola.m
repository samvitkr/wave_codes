close all
clear
baseDir = '/scratch.global/kuma0458/c8ak1_re180/run';
fign =                     'JApiola_c8ak1_re180.fig';
load(fullfile(baseDir, 'Piola_Relation8_LabFrame.mat'), 'times', 'dIa_dt', 'T_dt_a', 'T_Vn', 'T_R', 'RHS_total');
Lx=4*pi;
err =100*(Lx-RHS_total)./Lx

% Plot Diagnostic Balance
f=figure('Name', 'Piola JA Relation 8 Balance (Lab Frame)');
 subplot(1,2,1)
plot(times, dIa_dt, 'LineWidth', 1.5, 'DisplayName', 'dI_a/dt'); hold on;
plot(times, T_dt_a, 'LineWidth', 1.5, 'DisplayName', '\int u \cdot \partial_t a dV');
plot(times, T_Vn, 'LineWidth', 1.5, 'DisplayName', '\int V_n a \cdot u dS');
plot(times, T_R, 'LineWidth', 1.5, 'DisplayName', '\int a \cdot R dV');
plot(times, RHS_total, 'k-', 'LineWidth', 2, 'DisplayName', 'Total RHS ');
xlabel('Time');
ylabel('Virtual Work / Flux');
title('Piola Carrier Moving-Domain Diagnostic Balance');
legend('Location', 'best');
legend boxoff
yline(Lx,'DisplayName','K \Delta B')

grid on;

subplot(1,2,2)
plot(times,err,'-b')
ylabel("%error")
xlabel('Time')
saveas(f,fign)


