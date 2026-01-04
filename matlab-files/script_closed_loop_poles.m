%% @file script_closed_loop_poles.m
% @brief    Closed-Loop Stability and Pole-Zero Mapping for F4-E
%
% @details  This script evaluates the effectiveness of the robust design point Q1 
%           by mapping the closed-loop poles against the required stability 
%           regions (Gamma). It visualizes the MIL-F-8785C requirements as 
%           damping lines and frequency arcs (wa, wb), ensuring the aircraft 
%           short-period dynamics meet specific flying quality standards.

% --- SPECIFICATIONS ---
zeta_min = 0.35;
% Values from Military Specifications table
wa_vec = [2.02, 3.50, 2.19, 3.29]; 
wb_vec = [7.23, 12.6, 7.86, 11.8]; 

% --- VISUALIZATION ---
for i = 1:4
    figure('Name', ['Required Region - FC' num2str(i)]);
    hold on; grid on;

    % 1. Calculate Boundary Geometry
    phi = acos(zeta_min); 
    theta_range = linspace(pi - phi, pi + phi, 200);

    % Outer arc (omega_b)
    x_outer = wb_vec(i) * cos(theta_range);
    y_outer = wb_vec(i) * sin(theta_range);

    % Inner arc (omega_a)
    x_inner = wa_vec(i) * cos(fliplr(theta_range));
    y_inner = wa_vec(i) * sin(fliplr(theta_range));

    % 2. Create Shaded Region (Gamma_i)
    px = [x_outer, x_inner];
    py = [y_outer, y_inner];
    
    % Shading the region
    patch(px, py, [0.8 0.9 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');

    % 3. Plot Boundaries for Legend
    % Inner Arc
    plot(wa_vec(i) * cos(theta_range), wa_vec(i) * sin(theta_range), 'k:', 'LineWidth', 1.2, ...
        'DisplayName', ['$\omega_a = ', num2str(wa_vec(i)), '$ rad/s']);
    % Outer Arc
    plot(wb_vec(i) * cos(theta_range), wb_vec(i) * sin(theta_range), 'r--', 'LineWidth', 1.2, ...
        'DisplayName', ['$\omega_b = ', num2str(wb_vec(i)), '$ rad/s']);
    % Damping Lines
    w_ext = linspace(0, wb_vec(i)*1.5, 10);
    plot(-zeta_min*w_ext, w_ext*sqrt(1-zeta_min^2), 'b-.', 'LineWidth', 1, 'DisplayName', '$\zeta = 0.35$');
    plot(-zeta_min*w_ext, -w_ext*sqrt(1-zeta_min^2), 'b-.', 'LineWidth', 1, 'HandleVisibility', 'off');

    % 4. Center Label (Gamma_i)
    r_mid = (wa_vec(i) + wb_vec(i)) / 2;
    text(-r_mid, 0, ['$\Gamma_', num2str(i), '$'], 'Interpreter', 'latex', ...
        'FontSize', 18, 'HorizontalAlignment', 'center', 'Color', [0 0.2 0.5]);

    % 5. Apply User Preferred Plot Settings [cite: 2026-01-01]
    set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
    set(gca, 'XDir', 'normal', 'YDir', 'normal');
    
    xlabel('Real Axis', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Imaginary Axis', 'Interpreter', 'latex', 'FontSize', 12);
    title(['Required Pole Region: FC' num2str(i)], 'Interpreter', 'latex', 'FontSize', 14);
    
    % Axis and Legend
    axis equal;
    xlim([-wb_vec(i)*2 1]);
    ylim([-wb_vec(i)*1.1 wb_vec(i)*1.1]);
    legend('Location', 'best', 'Interpreter', 'latex');
end

%% --- DESIGN POINT & SPECIFICATIONS ---
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