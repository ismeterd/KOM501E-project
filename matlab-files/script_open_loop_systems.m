%% @file open_loop_system_analysis.m
% @brief    Open Loop System Responses of F4-E Aircraft
%
% @details  This script simulates short period longitudinal mode of an F4-E 
% aircraft with additional horizontal canards.
%
% @author   İsmet ERDOĞAN
% @date     2025-12-11

NUM_OF_FLIGHT_CONDS = 4;

% state system properties
stateNames = {  'normal accelaration'...        % n_z
                'pitch rate'...                 % q
                'elevator deflection'};         % delta_e
stateUnits = { ''...
               'rad/s'...
               'rad'};
inputName = "commanded elevator deflection";   % delta_e_comm
inputUnit = "rad";

% Flight Condition 1: 0.5  MACH -  5000' Altitude
% Flight Condition 2: 0.85 MACH -  5000' Altitude
% Flight Condition 3: 0.9  MACH - 35000' Altitude
% Flight Condition 4: 1.5  MACH - 35000' Altitude

% Flight Condition Table (FC1 | FC2 | FC3 | FC4)
a_11 = [-0.9896   -1.702    -0.667    -0.5162];
a_12 = [17.41     50.72     18.11     26.96];
a_13 = [96.15     263.5     84.34     178.9];

a_21 = [0.2648    0.2201    0.08201   -0.6896];
a_22 = [-0.8512   -1.418    -0.6587   -1.225];
a_23 = [-11.39    -31.99    -10.81    -30.38];

b_1  = [-97.78    -272.2    -85.09    -175.6];


% Preallocate struct array
FC(NUM_OF_FLIGHT_CONDS) = struct('A',[], 'B',[], 'C',[], 'D',[], ...
                                 'sys',[], 'name',[]);

for i=1:NUM_OF_FLIGHT_CONDS
    % create matrices and state space models
    A = [a_11(i), a_12(i), a_13(i);
         a_21(i), a_22(i), a_23(i);
               0,       0, -14];
    B = [b_1(i); 0; 14];
    C = eye(3);
    D = [0; 0; 0];

    fcName = sprintf("FC%d", i);
    nameString = "(Open Loop) F-4E at " + fcName;

    sys = ss(A, B, C, D, "Name", nameString);

    % system properties
    sys.StateName  = stateNames;
    sys.OutputName = stateNames;
    sys.InputName  = inputName;

    sys.StateUnit  = stateUnits;
    sys.OutputUnit = stateUnits;
    sys.InputUnit  = inputUnit;

    % Save into struct
    FC(i).A    = A;
    FC(i).B    = B;
    FC(i).C    = C;
    FC(i).D    = D;
    FC(i).sys  = sys;
    FC(i).name = nameString;
end

%% save data
for i = 1:NUM_OF_FLIGHT_CONDS
    fileName = sprintf("FC%d.mat", i);
    filePath = fullfile("data\", fileName);

    varName = sprintf("FC%d", i);
    S.(varName) = FC(i);

    save(filePath, "-struct", "S");

    clear S;
end