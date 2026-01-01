%% Ackermann Controller Parameter Space Analysis - All Flight Conditions
% This script computes and plots admissible gain regions for FC1, FC2, FC3, and FC4.

%% Section 1: Admissible Gain Region for Flight Condition 1
% Load FC1 data
fileName1 = "FC1.mat";
data1 = load(fullfile("data", fileName1));
A1 = data1.FC1.A; b1 = data1.FC1.B;
wa1 = 2.02; 

figure('Name', 'Admissible Gain Region for Flight Condition 1');
hold on; grid on;

% Boundary b-c: Damping Limit (zeta = 0.35)
w_vec = linspace(1, 40, 1500); 
knz_bc = []; kq_bc = [];

for w = w_vec
    s_val = -0.35*w + 1i*w*sqrt(1-0.35^2);
    [knz, kq] = solve_k_plane(s_val, A1, b1);
    
    if ~isnan(knz) && knz < 0 && kq < 0
        K_test = [knz, kq, 0]; 
        poles = eig(A1 - b1*K_test);
        real_poles = poles(abs(imag(poles)) < 1e-4);
        
        if any(real_poles <= -12.6 & real_poles >= -70.5)
            knz_bc = [knz_bc; knz]; 
            kq_bc = [kq_bc; kq];
        end
    end
end

plot(knz_bc, kq_bc, 'b', 'LineWidth', 2, 'DisplayName', 'b-c (\zeta=0.35)');

% Boundary a-b: Natural Frequency Limit (omega = 2.02)
zeta_vec = linspace(0.35, 1, 500);
knz_ab = []; kq_ab = [];

for z = zeta_vec
    s_val = -z*wa1 + 1i*wa1*sqrt(1-z^2);
    [knz, kq] = solve_k_plane(s_val, A1, b1);
    
    if ~isnan(knz) && knz < 0 && kq < 0
        knz_ab = [knz_ab; knz]; 
        kq_ab = [kq_ab; kq];
    end
end

plot(knz_ab, kq_ab, 'r', 'LineWidth', 2, 'DisplayName', ['a-b (\omega=',num2str(wa1),')']);

% Real Root Boundaries (sigma = -70 and -12.6)
knz_range = linspace(-0.4, 0.05, 1000);
sigmas = [-70, -12.6]; colors = {'g', 'm'};
names = {'c-d (\sigma=-70)', 'd-e (\sigma=-12.6)'};

for i = 1:length(sigmas)
    Phi = sigmas(i)*eye(3) - A1; v = adjoint(Phi)*b1;
    kq_line = -(det(Phi) + knz_range * v(1)) / v(2);
    mask = kq_line <= 0.1 & kq_line >= -3.5;
    
    if any(mask)
        plot(knz_range(mask), kq_line(mask), colors{i}, 'LineWidth', 1.5, 'DisplayName', names{i}); 
    end
end

% Plot Settings
xlim([-0.3 0]); ylim([-3 0.1]);
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'XDir', 'normal', 'YDir', 'normal');
xlabel('$k_{N_z}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$k_q$', 'Interpreter', 'latex', 'FontSize', 12);
title('Admissible Gain Region for Flight Condition 1', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'best');

%% Section 2: Admissible Gain Region for Flight Condition 2
% Load FC2 data
fileName2 = "FC2.mat";
data2 = load(fullfile("data", fileName2));
A2 = data2.FC2.A; b2 = data2.FC2.B;
wa2 = 3.50; 

figure('Name', 'Admissible Gain Region for Flight Condition 2');
hold on; grid on;

% Boundary b-c: Damping Limit (zeta = 0.35)
w_vec = linspace(1, 40, 1500); 
knz_bc = []; kq_bc = [];

for w = w_vec
    s_val = -0.35*w + 1i*w*sqrt(1-0.35^2);
    [knz, kq] = solve_k_plane(s_val, A2, b2);
    
    if ~isnan(knz) && knz < 0 && kq < 0
        K_test = [knz, kq, 0]; 
        poles = eig(A2 - b2*K_test);
        real_poles = poles(abs(imag(poles)) < 1e-4);
        
        if any(real_poles <= -12.6 & real_poles >= -70.5)
            knz_bc = [knz_bc; knz]; 
            kq_bc = [kq_bc; kq];
        end
    end
end

plot(knz_bc, kq_bc, 'b', 'LineWidth', 2, 'DisplayName', 'b-c (\zeta=0.35)');

% Boundary a-b: Natural Frequency Limit (omega = 3.50)
zeta_vec = linspace(0.35, 1, 500);
knz_ab = []; kq_ab = [];

for z = zeta_vec
    s_val = -z*wa2 + 1i*wa2*sqrt(1-z^2);
    [knz, kq] = solve_k_plane(s_val, A2, b2);
    
    if ~isnan(knz) && knz < 0 && kq < 0
        knz_ab = [knz_ab; knz]; 
        kq_ab = [kq_ab; kq];
    end
end

plot(knz_ab, kq_ab, 'r', 'LineWidth', 2, 'DisplayName', ['a-b (\omega=',num2str(wa2),')']);

% Real Root Boundaries
for i = 1:length(sigmas)
    Phi = sigmas(i)*eye(3) - A2; v = adjoint(Phi)*b2;
    kq_line = -(det(Phi) + knz_range * v(1)) / v(2);
    mask = kq_line <= 0.1 & kq_line >= -3.5;
    
    if any(mask)
        plot(knz_range(mask), kq_line(mask), colors{i}, 'LineWidth', 1.5, 'DisplayName', names{i}); 
    end
end

% Plot Settings
xlim([-0.3 0]); ylim([-3 0.1]);
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'XDir', 'normal', 'YDir', 'normal');
xlabel('$k_{N_z}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$k_q$', 'Interpreter', 'latex', 'FontSize', 12);
title('Admissible Gain Region for Flight Condition 2', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'best');

%% Section 3: Admissible Gain Region for Flight Condition 3
% Load FC3 data
fileName3 = "FC3.mat";
data3 = load(fullfile("data", fileName3));
A3 = data3.FC3.A; b3 = data3.FC3.B;
wa3 = 2.19; 

figure('Name', 'Admissible Gain Region for Flight Condition 3');
hold on; grid on;

% Boundary b-c: Damping Limit (zeta = 0.35)
w_vec = linspace(1, 40, 1500); 
knz_bc = []; kq_bc = [];

for w = w_vec
    s_val = -0.35*w + 1i*w*sqrt(1-0.35^2);
    [knz, kq] = solve_k_plane(s_val, A3, b3);
    
    if ~isnan(knz) && knz < 0 && kq < 0
        K_test = [knz, kq, 0]; 
        poles = eig(A3 - b3*K_test);
        real_poles = poles(abs(imag(poles)) < 1e-4);
        
        if any(real_poles <= -12.6 & real_poles >= -70.5)
            knz_bc = [knz_bc; knz]; 
            kq_bc = [kq_bc; kq];
        end
    end
end

plot(knz_bc, kq_bc, 'b', 'LineWidth', 2, 'DisplayName', 'b-c (\zeta=0.35)');

% Boundary a-b: Natural Frequency Limit (omega = 2.19)
zeta_vec = linspace(0.35, 1, 500);
knz_ab = []; kq_ab = [];

for z = zeta_vec
    s_val = -z*wa3 + 1i*wa3*sqrt(1-z^2);
    [knz, kq] = solve_k_plane(s_val, A3, b3);
    
    if ~isnan(knz) && knz < 0 && kq < 0
        knz_ab = [knz_ab; knz]; 
        kq_ab = [kq_ab; kq];
    end
end

plot(knz_ab, kq_ab, 'r', 'LineWidth', 2, 'DisplayName', ['a-b (\omega=',num2str(wa3),')']);

% Real Root Boundaries
for i = 1:length(sigmas)
    Phi = sigmas(i)*eye(3) - A3; v = adjoint(Phi)*b3;
    kq_line = -(det(Phi) + knz_range * v(1)) / v(2);
    mask = kq_line <= 0.1 & kq_line >= -3.5;
    
    if any(mask)
        plot(knz_range(mask), kq_line(mask), colors{i}, 'LineWidth', 1.5, 'DisplayName', names{i}); 
    end
end

% Plot Settings
xlim([-0.3 0]); ylim([-3 0.1]);
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'XDir', 'normal', 'YDir', 'normal');
xlabel('$k_{N_z}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$k_q$', 'Interpreter', 'latex', 'FontSize', 12);
title('Admissible Gain Region for Flight Condition 3', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'best');

%% Section 4: Admissible Gain Region for Flight Condition 4
% Load FC4 data
fileName4 = "FC4.mat";
data4 = load(fullfile("data", fileName4));
A4 = data4.FC4.A; b4 = data4.FC4.B;
wa4 = 3.29; 

figure('Name', 'Admissible Gain Region for Flight Condition 4');
hold on; grid on;

% Boundary b-c: Damping Limit (zeta = 0.35)
w_vec = linspace(1, 120, 5000); 
knz_bc = []; kq_bc = [];

for w = w_vec
    s_val = -0.35*w + 1i*w*sqrt(1-0.35^2);
    [knz, kq] = solve_k_plane(s_val, A4, b4);
    
    if ~isnan(knz) && knz < 0.05 && kq < 0.05
        K_test = [knz, kq, 0]; 
        poles = eig(A4 - b4*K_test);
        real_poles = poles(abs(imag(poles)) < 1e-4);
        
        if any(real_poles <= -12.5 & real_poles >= -70.5)
            knz_bc = [knz_bc; knz]; 
            kq_bc = [kq_bc; kq];
        end
    end
end

plot(knz_bc, kq_bc, 'b', 'LineWidth', 2, 'DisplayName', 'b-c (\zeta=0.35)');

% Boundary a-b: Natural Frequency Limit (omega = 3.29)
zeta_vec = linspace(0.1, 1, 1500); 
knz_ab = []; kq_ab = [];

for z = zeta_vec
    s_val = -z*wa4 + 1i*wa4*sqrt(1-z^2);
    [knz, kq] = solve_k_plane(s_val, A4, b4);
    
    if ~isnan(knz) && knz <= 0.05 && kq <= 0.05
        knz_ab = [knz_ab; knz]; 
        kq_ab = [kq_ab; kq];
    end
end

if ~isempty(knz_ab), plot(knz_ab, kq_ab, 'r', 'LineWidth', 2, 'DisplayName', ['a-b (\omega=',num2str(wa4),')']); end

% Real Root Boundaries
knz_range4 = linspace(-0.6, 0.1, 2500);

for i = 1:length(sigmas)
    Phi = sigmas(i)*eye(3) - A4; v = adjoint(Phi)*b4;
    kq_line = -(det(Phi) + knz_range4 * v(1)) / v(2);
    mask = kq_line <= 0.1 & kq_line >= -3.5;
    
    if any(mask)
        plot(knz_range4(mask), kq_line(mask), colors{i}, 'LineWidth', 1.5, 'DisplayName', names{i}); 
    end
end

% Plot Settings
xlim([-0.4 0]); ylim([-3.5 0.1]);
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'XDir', 'normal', 'YDir', 'normal');
xlabel('$k_{N_z}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$k_q$', 'Interpreter', 'latex', 'FontSize', 12);
title('Admissible Gain Region for Flight Condition 4', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'best');

%% Section 5: Robust Admissible Region (Rnom) - Final Intersection
% Dynamically load all .mat files into a single struct array
for k = 1:4
    temp_data = load(fullfile("data", sprintf("FC%d.mat", k)));
    FC(k) = temp_data.(sprintf("FC%d", k));
end

figure('Name', 'Robust Admissible Region - Intersection');
hold on; grid on;

% 1. Plot Damping Boundaries (zeta = 0.35) for all 4 Flight Conditions
colors = {'r', 'b', 'g', 'm'};

for i = 1:4
    A = FC(i).A; b = FC(i).B;
    w_vec = linspace(1, 60, 3000); 
    knz_vals = []; kq_vals = [];
    
    for w = w_vec
        s_val = -0.35*w + 1i*w*sqrt(1-0.35^2);
        [knz, kq] = solve_k_plane(s_val, A, b);
        
        if ~isnan(knz) && knz < 0.1 && kq < 0.1
            K_test = [knz, kq, 0]; 
            poles = eig(A - b*K_test);
            % Ghost line filter: Ensure the real pole is within stability limits
            real_poles = poles(abs(imag(poles)) < 1e-4);
            
            if any(real_poles <= -12.5 & real_poles >= -70.5)
                knz_vals = [knz_vals; knz]; 
                kq_vals = [kq_vals; kq];
            end
        end
    end
    
    plot(knz_vals, kq_vals, colors{i}, 'LineWidth', 1.5, 'DisplayName', ['FC' num2str(i) ' \zeta=0.35']);
end

% 2. Plot Critical Sigma Boundaries as SOLID Lines
knz_range = linspace(-0.4, 0.1, 1000);

% a) FC2 Actuator Limit: sigma = -70 (Solid Black)
v2 = adjoint(-70*eye(3) - FC(2).A) * FC(2).B;
kq_s70 = -(det(-70*eye(3)-FC(2).A) + knz_range * v2(1)) / v2(2);
plot(knz_range, kq_s70, 'k', 'LineWidth', 2, 'DisplayName', '\sigma_2 = -70');

% b) FC1 Stability Limit: sigma = -7.23 (Solid Dark Gray)
v1_low = adjoint(-7.23*eye(3) - FC(1).A) * FC(1).B;
kq_s723 = -(det(-7.23*eye(3)-FC(1).A) + knz_range * v1_low(1)) / v1_low(2);
plot(knz_range, kq_s723, 'Color', [0.3 0.3 0.3], 'LineWidth', 2, 'DisplayName', '\sigma_1 = -7.23');

% c) FC1 Speed/Landing Limit: sigma = -2.02 (Solid Light Gray)
v1_high = adjoint(-2.02*eye(3) - FC(1).A) * FC(1).B;
kq_s202 = -(det(-2.02*eye(3)-FC(1).A) + knz_range * v1_high(1)) / v1_high(2);
plot(knz_range, kq_s202, 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'DisplayName', '\sigma_1 = -2.02');

% 3. Mark Design Point Q1
plot(-0.115, -0.8, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'Design Point Q_1');
text(-0.11, -0.7, '$Q_1$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');

% Plot Settings
xlim([-0.25 0]); ylim([-2 0]);
set(gca, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'XDir', 'normal', 'YDir', 'normal');
xlabel('$k_{N_z}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$k_q$', 'Interpreter', 'latex', 'FontSize', 12);
title('Robust Admissible Region for All Flight Conditions ($R_{nom}$)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Location', 'best');

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