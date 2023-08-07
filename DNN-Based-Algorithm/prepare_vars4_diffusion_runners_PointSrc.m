function [ tx_node, rx_node, env_params, sim_params ] = prepare_vars4_diffusion_runners_PointSrc(dist, emission_pt, receiver_pt, r_r, D, D_tx, D_rx, delta_t, molecules_perTs, ts_inSeconds, tss_inSeconds, symbol_probs, nsym, replication)
% Description:
% set the simulator parameters
dist_inMicroMeters = dist; % distance
tx_node.emission_point     = emission_pt; % the coordinate of transmitter
tx_node.D_inMicroMeterSqrPerSecond = D_tx; % µm^2/s
rx_node.r_inMicroMeters    = r_r; % µm
rx_node.center             = receiver_pt; % the coordinate of receiver
rx_node.synch_offset       = 0; % Param
rx_node.p_react            = 4; % p_react = 1 perfect_absorption p_react = 2 imperfect_absorption p_react = 3 mobile perfect_absorption
rx_node.D_inMicroMeterSqrPerSecond = D_rx; % the diffusion coefficient of receiver
env_params.snr_db                         = 30; % NOISE FROM Environment  
env_params.D_inMicroMeterSqrPerSecond     = D; % the diffusion coefficient of molecules
env_params.destruction_limit              = 7*dist_inMicroMeters + 45 + rx_node.r_inMicroMeters; % the maximum distance

sim_params.delta_t                        = delta_t; % individual sample time in the simulation time
sim_params.molecules_perTs                = molecules_perTs;   % Param changes in Runner
sim_params.ts_inSeconds                   = ts_inSeconds; % Param changes in Runner
sim_params.tss_inSeconds                  = tss_inSeconds; % Param changes in Runner

sim_params.symbol_probs    = symbol_probs;
sim_params.nsym            = nsym; 
sim_params.replication     = replication; % Increase Replication


end
