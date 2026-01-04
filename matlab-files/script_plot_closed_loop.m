%% @file script_plot_closed_loop.m
% @brief    Closed-Loop Simulation and Performance Analysis of F4-E
%
% @details  This script executes Simulink models for all four flight conditions (FC1-FC4)
%           using the robust design point Q1. It plots c* step responses.

mdl = "simModel_closed_loop_systems";
blk = mdl + "/F-4E with Canards";

choices = ["FC1","FC2","FC3","FC4"];

for c = choices
    in = Simulink.SimulationInput(mdl);
    in = in.setBlockParameter(blk, "FCSelection", c);
    in = in.setModelParameter("StopTime", "2");
    simOut = sim(in);

    % post-process...
    output = simOut.output.Data;
    delta = simOut.delta.Data;
    t = simOut.tout;
    
    % Process the simulation output for further analysis
    figure;
    plot(t, output./output(end), "LineWidth", 2);
    xlabel('time $(s)$', "Interpreter", "latex");
    ylabel('$c^*(t)$', "Interpreter", "latex");
    title(sprintf("Closed-Loop Response at %s", c), ...
        "Interpreter", "latex", "FontSize", 14);
    grid on;
end

