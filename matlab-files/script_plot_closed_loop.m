mdl = "simModel_closed_loop_systems";
blk = mdl + "/F-4E with Canards";

choices = ["FC1","FC2","FC3","FC4"];

for c = choices
    in = Simulink.SimulationInput(mdl);
    in = in.setBlockParameter(blk, "FCSelection", c);
    simOut = sim(in);

    % post-process...
    output = simOut.output.Data;
    delta = simOut.delta.Data;
    t = simOut.tout;
    
    % Process the simulation output for further analysis
    figure;
    plot(t, output./output(end));
    xlabel('Time (s)');
    ylabel('Output Response');
    title('Closed Loop System Response');
    grid on;
end

