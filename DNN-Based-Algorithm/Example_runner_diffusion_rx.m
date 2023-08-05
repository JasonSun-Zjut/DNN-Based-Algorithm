function [ b_Rx_1, b_Rx_2] = Example_runner_diffusion_rx(n, NA,coordinate_vector,r,diffusion_coefficient)
% Description:
% NA : number of molecules to emit
% coordinate_vector: [[Tx1];[Rx1];[Rx2]]
% r : radius of Rx
% diffusion_coefficient : [D_mol,D_tx,D_rx]
% n: nth time slot emit molecules
%% experiment
fprintf(1, '\n################# 开始模拟Tx在 %dth时隙释放分子的过程  ############', n);
delta_t              = 0.001;
fprintf(1,'\n ## delta_t = %f s', delta_t);
num_molecules_to_emit = NA;
fprintf(1,'\n ## Num. Emitted Molecules = %d ', num_molecules_to_emit);
nsym                 = 1;
ts                   = 5;
% Prepare Variables
emission_pt = coordinate_vector(1,:);
receiver_pt = coordinate_vector(2,:);
dist_inMicroMeters = norm(receiver_pt(1,:)-emission_pt(1,:))-r;
%dist_inMicroMeters, emission_pt, receiver_pt, r_r, D, D_tx, D_rx, delta_t, molecules_perTs, ts_inSeconds, tss_inSeconds, symbol_probs, nsym, replication%
[tx_node, rx_node, env_params, sim_params] = prepare_vars4_diffusion_runners_PointSrc(dist_inMicroMeters, emission_pt, receiver_pt, r, diffusion_coefficient(1), diffusion_coefficient(2), diffusion_coefficient(3), delta_t, num_molecules_to_emit, ts, 0.001, [0.5 0.5], nsym, 200);
% Run 
res_Signal = runner_diffusion_rx(n, tx_node, rx_node, env_params, sim_params);
end

function [res] = runner_diffusion_rx(n, tx_node, rx_node, env_params, sim_params)
tx_sym_matrix              = repmat([zeros(1,n-1),1], sim_params.replication, 1);

tx_node.mod                = 0; %% BCSK (pulse)
rx_node.demod              = tx_node.mod;

[ nRx_raw_matrix_wout_noise, stats ] = CORE_sim_replicator( ...
           'CORE_sim_diffusion_3d_P2S_wAbsorption', ...% Simulator Name ##You Can Change This##
           tx_sym_matrix, ... % Tx Symbol Sequences
           tx_node,    ...% Tx node properties
           rx_node,    ...% Rx node properties
           env_params, ...% Environment properties
           sim_params );

res.nRx_raw_matrix_wout_noise = nRx_raw_matrix_wout_noise;
res.stats = stats;
ts_step       =  round( sim_params.ts_inSeconds / sim_params.delta_t );
res.nRx_1_avg =  sum(nRx_raw_matrix_wout_noise(:,(n-1)*ts_step+1:n*ts_step,:), 3) / size(nRx_raw_matrix_wout_noise(:,(n-1)*ts_step+1:n*ts_step,:), 3);
%res.nRx_2_avg =  sum(nRx_raw_matrix_wout_noise(:,(end-ts_step)+1:end,:), 3) / size(nRx_raw_matrix_wout_noise(:,(end-ts_step)+1:end,:), 3);
end