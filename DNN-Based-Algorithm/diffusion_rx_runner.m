function [res_Signal, b_Rx_1, b_Rx_2, b_Rx_3] = diffusion_rx_runner(NA, coordinate, diffusion_coefficient)
% Description:
% the file to run molecular diffusion simulator
% NA : number of molecules to emit
% coordinate_vector: [[Tx1];[Tx2];[Tx2];[Rx1]]
% r : radius of Rx
% diffusion_coefficient : [D_mol,D_tx,D_rx]
% n: nth time slot emit molecules

% experiment
r = 3; 
n = 5;
fprintf(1, '\n################# 开始模拟Tx在 %dth时隙释放分子的过程  ############', n);
delta_t              = 0.001;
fprintf(1,'\n ## delta_t = %f s', delta_t);
num_molecules_to_emit = NA;
fprintf(1,'\n ## Num. Emitted Molecules = %d ', num_molecules_to_emit);
nsym                 = 1;
ts                   = 5;
% Prepare Variables
transmitter_pt = coordinate(1:3,:);
receiver_pt = coordinate(4,:);
dist_inMicroMeters = norm(repmat(receiver_pt(1,:), 1, 3)-transmitter_pt(1:3,:)) - r;
[tx_node, rx_node, env_params, sim_params] = prepare_vars4_diffusion_runners_PointSrc(dist_inMicroMeters, transmitter_pt, receiver_pt, r, diffusion_coefficient(1), diffusion_coefficient(2), diffusion_coefficient(3), delta_t, num_molecules_to_emit, ts, 0.001, [0.5 0.5], nsym, 200);
% Run 
res_Signal = run_diffusion_rx(n, tx_node, rx_node, env_params, sim_params);
signal_resolution_merge = 1;
x_time = 1e-10:(delta_t*signal_resolution_merge):ts;
% non_linear_model
nRx_1 = res_Signal.nRx_1_avg;
nRx_2 = res_Signal.nRx_2_avg;
nRx_3 = res_Signal.nRx_3_avg;
y_1 = cumsum(nRx_1)./NA;
y_2 = cumsum(nRx_2)./NA;
y_3 = cumsum(nRx_3)./NA;
V_rx = 4*pi*r^3/3;
d_tx_rx_1 = norm(receiver_pt(1,:)-emission_pt)-r;
d_tx_rx_2 = norm(receiver_pt(2,:)-emission_pt)-r;
D1 = diffusion_coefficient(1) + diffusion_coefficient(3);
D2 = diffusion_coefficient(2) + diffusion_coefficient(3);
model_f_Rx_1 = @(b,x)b(1)*V_rx.*erfc(d_tx_rx_1./((4*D2)^b(2).*(D1/D2.*x+(n-1)*ts).^b(3)))./(4*pi*D2*d_tx_rx_1);
model_f_Rx_2 = @(b,x)b(1)*V_rx.*erfc(d_tx_rx_2./((4*D2)^b(2).*(D1/D2.*x+(n-1)*ts).^b(3)))./(4*pi*D2*d_tx_rx_2);
model_f_Rx_3 = @(b,x)b(1)*V_rx.*erfc(d_tx_rx_3./((4*D2)^b(2).*(D1/D2.*x+(n-1)*ts).^b(3)))./(4*pi*D2*d_tx_rx_3);
b_init = 0.5 + (1 - 0.5).*rand(3,1);
Rx_1_mdl = fitnlm(x_time(1:end),y_1(1:end),model_f_Rx_1,b_init);
Rx_2_mdl = fitnlm(x_time(1:end),y_2(1:end),model_f_Rx_2,b_init);
Rx_3_mdl = fitnlm(x_time(1:end),y_3(1:end),model_f_Rx_3,b_init);
% fitted parameter b1, b2, b3
b_Rx_1 = round(Rx_1_mdl.Coefficients.Estimate, 3);
b_Rx_2 = round(Rx_2_mdl.Coefficients.Estimate, 3);
b_Rx_3 = round(Rx_3_mdl.Coefficients.Estimate, 3);
end

function [res] = run_diffusion_rx(n, tx_node, rx_node, env_params, sim_params)
tx_sym_matrix              = repmat([zeros(1,n-1),1], sim_params.replication, 1);

tx_node.mod                = 0; %% BCSK (pulse)
rx_node.demod              = tx_node.mod;

[ nRx_raw_matrix_wout_noise, stats ] = CORE_sim_replicator( ...
           'CORE_sim_diffusion_3d_P2S_wAbsorption', ...% Simulator Function
           tx_sym_matrix, ... % Tx Symbol Sequences
           tx_node,    ...% Tx node properties
           rx_node,    ...% Rx node properties
           env_params, ...% Environment properties
           sim_params );

res.nRx_raw_matrix_wout_noise = nRx_raw_matrix_wout_noise;
res.stats = stats;
ts_step       =  round( sim_params.ts_inSeconds / sim_params.delta_t ); % the number of samples
res.nRx_1_avg =  sum(nRx_raw_matrix_wout_noise(:,(n-1) * ts_step + 1:n * ts_step,:), 3) / size(nRx_raw_matrix_wout_noise(:,(n-1)*ts_step+1:n*ts_step,:), 3); % the number of molecules Rx1 observed in each sample time
res.nRx_2_avg =  sum(nRx_raw_matrix_wout_noise(:,((end - 2) * ts_step) + 1:(end - 1) * ts_step,:), 3) / size(nRx_raw_matrix_wout_noise(:,(end - 2) * ts_step + 1 :(end-1) * ts_step,:), 3); % the number of molecules Rx2 observed in each sample time
res.nRx_3_avg =  sum(nRx_raw_matrix_wout_noise(:,(end - ts_step) + 1:end,:), 3) / size(nRx_raw_matrix_wout_noise(:,(end - ts_step) + 1:end,:), 3); % the number of molecules Rx3 observed in each sample time
end