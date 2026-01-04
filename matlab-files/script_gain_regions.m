%% @file script_gain_regions.m
% @brief    Robust Admissible Gain Region Analysis for F4-E Aircraft
%
% @details  This script performs a parameter space analysis using the D-Decomposition 
%           method to determine the robust admissible controller gains (knz, kq). 
%           It incorporates MIL-F-8785C flying quality requirements including 
%           damping ratios, natural frequency boundaries, and structural limits 
%           across four critical flight conditions (FC1-FC4).

sigma_actuator = -70; % Structural/Actuator frequency limit

% Flight Data Matrix: [wa, wb, FC_index]
fc_params = [2.02, 7.23, 1;   % FC1
             3.50, 12.6, 2;   % FC2
             2.19, 7.86, 3;   % FC3
             3.29, 11.8, 4];  % FC4

for i = 1:size(fc_params, 1)
    wa = fc_params(i,1); wb = fc_params(i,2); idx = fc_params(i,3);
    
    % Load flight condition data
    data = load(fullfile("data", sprintf("FC%d.mat", idx)));
    A = data.(sprintf("FC%d", idx)).A; b = data.(sprintf("FC%d", idx)).B;
    
    figure('Name', ['Admissible Gain Region for FC', num2str(idx)]);
    hold on; grid on;
    
    % --- 1. Boundary: Damping Limit (zeta = 0.35) with Ghost Line Filter ---
    % Filtering ensures the third (real) pole stays within MIL-spec limits
    w_vec = linspace(wa, wb, 3000); 
    knz_bc = []; kq_bc = [];
    for w = w_vec
        s_val = -0.35*w + 1i*w*sqrt(1-0.35^2);
        [knz, kq] = solve_k_plane(s_val, A, b);
        
        if ~isnan(knz) && knz < 0.1 && kq < 0.1
            K_test = [knz, kq, 0]; 
            poles = eig(A - b*K_test);
            real_poles = poles(abs(imag(poles)) < 1e-4);
            
            % Check if the parasitic real pole is between -wb and -70
            if any(real_poles <= -wb & real_poles >= sigma_actuator)
                knz_bc = [knz_bc; knz]; 
                kq_bc = [kq_bc; kq];
            end
        end
    end
    plot(knz_bc, kq_bc, 'b', 'LineWidth', 2, 'DisplayName', '\zeta=0.35');

    % --- 2. Boundary: Lower Natural Frequency Arc (omega_a) ---
    zeta_vec = linspace(0.35, 1, 1000);
    knz_wa = []; kq_wa = [];
    for z = zeta_vec
        s_val = -z*wa + 1i*wa*sqrt(1-z^2);
        [knz, kq] = solve_k_plane(s_val, A, b);
        if ~isnan(knz), knz_wa = [knz_wa; knz]; kq_wa = [kq_wa; kq]; end
    end
    plot(knz_wa, kq_wa, 'r', 'LineWidth', 2, 'DisplayName', ['\omega_a=', num2str(wa)]);

    % --- 3. Boundary: Upper Natural Frequency Arc (omega_b) ---
    knz_wb_arc = []; kq_wb_arc = [];
    for z = zeta_vec
        s_val = -z*wb + 1i*wb*sqrt(1-z^2);
        [knz, kq] = solve_k_plane(s_val, A, b);
        if ~isnan(knz), knz_wb_arc = [knz_wb_arc; knz]; kq_wb_arc = [kq_wb_arc; kq]; end
    end
    plot(knz_wb_arc, kq_wb_arc, 'm', 'LineWidth', 2, 'DisplayName', ['\omega_b=', num2str(wb)]);

    % --- 4. Boundary: Real Root Boundaries (RRB) ---
    knz_range = linspace(-0.5, 0.1, 2000);
    % Sigma = -wb and Sigma = -70 (Structural limit)
    sigmas = [-wb, sigma_actuator];
    sigma_colors = {'m', 'g'};
    sigma_labels = {['\sigma=-\omega_b'], '\sigma=-70'};
    
    for s_idx = 1:length(sigmas)
        Phi = sigmas(s_idx)*eye(3) - A; v = adjoint(Phi)*b;
        kq_line = -(det(Phi) + knz_range * v(1)) / v(2);
        mask = kq_line <= 0.1 & kq_line >= -4;
        if any(mask), plot(knz_range(mask), kq_line(mask), 'Color', sigma_colors{s_idx}, ...
            'LineStyle', '--', 'LineWidth', 1.2, 'DisplayName', sigma_labels{s_idx}); end
    end

    % Applying plot settings per user preference
    set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'XDir', 'normal', 'YDir', 'normal');
    xlabel('$k_{N_z}$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$k_q$', 'Interpreter', 'latex', 'FontSize', 12);
    title(['Admissible Gain Region for FC', num2str(idx)], 'Interpreter', 'latex', 'FontSize', 14);
    xlim([-0.4 0.05]); ylim([-3.5 0.1]); legend('Location', 'best');
end

%% Section 5: Robust Admissible Region (Rnom) - Final Visualization with Values
% This script displays the intersection of all FC boundaries with frequency values in legend.

figure('Name', 'Robust Admissible Region - All FC Boundaries');
hold on; grid on;

% 1. Universal and Adaptive Constraints
sigma_actuator = -70; 
wa_vec = [2.02, 3.50, 2.19, 3.29]; % sigma_min limits for FC1-4
wb_vec = [7.23, 12.6, 7.86, 11.8]; % sigma_max limits for FC1-4
fc_colors = {'r', 'b', 'g', 'm'};    

% Ensure FC data structure is loaded
for k = 1:4
    temp = load(fullfile("data", sprintf("FC%d.mat", k)));
    FC_int(k) = temp.(sprintf("FC%d", k));
end

% 2. Plotting loop for all 4 Flight Conditions
knz_range = linspace(-0.3, 0.05, 1500);

for i = 1:4
    A = FC_int(i).A; b = FC_int(i).B;
    curr_c = fc_colors{i};
    wa = wa_vec(i); wb = wb_vec(i);
    
    % --- A) Damping Boundary (zeta = 0.35) ---
    knz_z = []; kq_z = [];
    for w = linspace(wa, wb, 1500)
        s_val = -0.35*w + 1i*w*sqrt(1-0.35^2);
        [knz, kq] = solve_k_plane(s_val, A, b);
        if ~isnan(knz) && knz < 0.1 && kq < 0.1
            % Ghost-line filter check
            K_t = [knz, kq, 0]; poles = eig(A - b*K_t);
            real_p = poles(abs(imag(poles)) < 1e-4);
            if any(real_p <= -wb & real_p >= sigma_actuator)
                knz_z = [knz_z; knz]; kq_z = [kq_z; kq];
            end
        end
    end
    % Displaying zeta and frequency range in legend
    plot(knz_z, kq_z, 'Color', curr_c, 'LineWidth', 2, ...
        'DisplayName', ['FC', num2str(i), ': \zeta=0.35, \omega \in [', num2str(wa), ',', num2str(wb), ']']);

    % --- B) Sigma_b Boundary (sigma = -wb) ---
    v_b = adjoint(-wb*eye(3) - A) * b;
    kq_wb = -(det(-wb*eye(3)-A) + knz_range * v_b(1)) / v_b(2);
    plot(knz_range, kq_wb, 'Color', curr_c, 'LineStyle', '--', 'LineWidth', 1.2, ...
        'DisplayName', ['FC', num2str(i), ': \sigma_b = -', num2str(wb)]);

    % --- C) Sigma_a Boundary (sigma = -wa) ---
    v_a = adjoint(-wa*eye(3) - A) * b;
    kq_wa = -(det(-wa*eye(3)-A) + knz_range * v_a(1)) / v_a(2);
    plot(knz_range, kq_wa, 'Color', curr_c, 'LineStyle', ':', 'LineWidth', 1.2, ...
        'DisplayName', ['FC', num2str(i), ': \sigma_a = -', num2str(wa)]);
end

% 3. Global Structural Limit (sigma = -70)
v_act = adjoint(sigma_actuator*eye(3) - FC_int(2).A) * FC_int(2).B;
kq_s70 = -(det(sigma_actuator*eye(3)-FC_int(2).A) + knz_range * v_act(1)) / v_act(2);
plot(knz_range, kq_s70, 'k', 'LineWidth', 2.5, 'DisplayName', '\sigma_{actuator} = -70');

% 4. Mark Robust Design Point Q1
plot(-0.115, -0.8, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'DisplayName', 'Design Point $Q_1$');
text(-0.11, -0.7, '$Q_1$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');

% Final Plot Settings [User Custom Preference]
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'XDir', 'normal', 'YDir', 'normal');
xlabel('$k_{N_z}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$k_q$', 'Interpreter', 'latex', 'FontSize', 12);
title('Robust Admissible Region ($R_{nom}$): Comprehensive Analysis', 'Interpreter', 'latex', 'FontSize', 14);
xlim([-0.25 0]); ylim([-1.75 0]); legend('Location', 'bestoutside', 'FontSize', 8);

%% Helper Function
function [knz, kq] = solve_k_plane(s, A, b)
    Phi = s*eye(3) - A; v = adjoint(Phi)*b;
    M = [real(v(1)), real(v(2)); imag(v(1)), imag(v(2))];
    Y = [-real(det(Phi)); -imag(det(Phi))];
    if abs(det(M)) > 1e-12
        sol = M \ Y; knz = sol(1); kq = sol(2);
    else
        knz = NaN; kq = NaN;
    end
end