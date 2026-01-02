%% Closed-Loop Pole Analysis with MIL-Spec Frequency Arcs
% This script visualizes the CL poles for each FC at Q1, showing 
% damping and natural frequency boundaries as per military standards.

% --- DESIGN POINT & SPECIFICATIONS ---
knz_q1 = -0.115;
kq_q1 = -0.8;
k_delta = 0;
K_q1 = [knz_q1, kq_q1, k_delta]; 

% Frequency limits from Military Specifications
wa_vec = [2.02, 3.50, 2.19, 3.29]; 
wb_vec = [7.23, 12.6, 7.86, 11.8]; 
zeta_min = 0.35;

% --- DATA PREPARATION ---
for k = 1:4
    temp_data = load(fullfile("data", sprintf("FC%d.mat", k)));
    FC(k) = temp_data.(sprintf("FC%d", k));
end

% --- INDIVIDUAL POLE-ZERO MAPS ---
for i = 1:4
    figure('Name', ['Closed-Loop Analysis - FC' num2str(i)]);
    hold on; grid on;
    
    % 1. Define Closed-Loop System and Poles
    A_cl = FC(i).A - FC(i).B * K_q1;
    cl_poles = eig(A_cl);
    
    % 2. Draw Military Specification Boundaries (Arcs and Lines)
    theta = linspace(pi/2, 3*pi/2, 500); % Arc for the left half-plane
    
    % a) Lower Frequency Arc (wa)
    plot(wa_vec(i)*cos(theta), wa_vec(i)*sin(theta), 'k:', 'LineWidth', 1.2, ...
        'DisplayName', ['\omega_a = ', num2str(wa_vec(i))]);
    
    % b) Upper Frequency Arc (wb)
    plot(wb_vec(i)*cos(theta), wb_vec(i)*sin(theta), 'r--', 'LineWidth', 1.2, ...
        'DisplayName', ['\omega_b = ', num2str(wb_vec(i))]);
    
    % c) Damping Ratio Lines (zeta = 0.35)
    w_line = linspace(0, 40, 500);
    plot(-zeta_min*w_line, w_line*sqrt(1-zeta_min^2), 'b:', 'LineWidth', 1.2, ...
        'DisplayName', '\zeta = 0.35');
    plot(-zeta_min*w_line, -w_line*sqrt(1-zeta_min^2), 'b:', 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    
    % 3. Plot Closed-Loop Poles Manually (To fix Legend and Marker Size)
    % Using large 'x' markers as requested
    plot(real(cl_poles), imag(cl_poles), 'bx', 'MarkerSize', 8, 'LineWidth', 2.5, ...
        'DisplayName', 'CL Poles');
    
    % 4. Formatting & Axis Settings (Interpreter Fixed)
    xlabel('Real Axis ($sec^{-1}$)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Imaginary Axis ($rad/sec$)', 'Interpreter', 'latex', 'FontSize', 12);
    title(['Closed-Loop Pole Locations of FC' num2str(i)], 'Interpreter', 'latex', 'FontSize', 14);
    
    axis equal;
    xlim([-40 0]); 
    ylim([-15 15]);
    
    legend('Location', 'best');
end