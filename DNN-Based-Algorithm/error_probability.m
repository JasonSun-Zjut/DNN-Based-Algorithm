function [Pe] = error_probability(n, NA, distance, r, diffusion_coefficient)
% Description: error probability of Rx in current time slot
% n: nth time slot
% NA: number of molecules to emit

emit_point = [0, 0];
rx_1_point = rand_coordinate_generate(emit_point, distance(1));
rx_2_point = rand_coordinate_generate(emit_point, distance(2));
coordinate_vector = [[emit_point,0];[rx_1_point,0];[rx_2_point,0]];
fprintf(1, '\n############ TX_1:[%d, %d, %d],Rx_1:[%d, %d, %d],Rx_2:[%d, %d, %d] ############\n',coordinate_vector(1,1:end),coordinate_vector(2,1:end),coordinate_vector(3,1:end));
%% experiment and non_linear_model
[Rx1_estimate_coefficient, Rx2_estimate_coefficient] = Example_runner_diffusion_rx(1, NA, coordinate_vector, r, diffusion_coefficient);

%% numberial analytical
Tx_1_coordinate = coordinate_vector(1,:);
Rx_1_coordinate = coordinate_vector(2,:);
Rx_2_coordinate = coordinate_vector(3,:);  
distance = [norm(Rx_1_coordinate-Tx_1_coordinate)-r,norm(Rx_2_coordinate-Tx_1_coordinate)-r];
T = 1;
% sample time 
tau = 0.005;

%% Tx_1 transmit bit 0
mu_0 = 0;
for k = 1:1:n-1
[Rx1_estimate_coefficient, Rx2_estimate_coefficient] = Example_runner_diffusion_rx(k, NA, coordinate_vector, r, diffusion_coefficient);
probability_function = probability_mobile(Rx1_estimate_coefficient, Rx2_estimate_coefficient, r, diffusion_coefficient, (k-1)*T,(n-k)*T+tau,distance);
fprintf(1, '\n############ mobiled TX_1----RX_1 probability ############\n');
disp(probability_function(1));
fprintf(1, '\n############ mobiled TX_1----RX_2 probability ############\n');
disp(probability_function(2));
mu_0 = mu_0 + 0.5 * NA * probability_function(1) +  0.5 * NA * probability_function(2);
end
[Rx1_estimate_coefficient, Rx2_estimate_coefficient] = Example_runner_diffusion_rx(n, NA, coordinate_vector, r, diffusion_coefficient);
probability_function = probability_mobile(Rx1_estimate_coefficient, Rx2_estimate_coefficient, r, diffusion_coefficient, (n-1)*T, tau, distance);
mu_0 = 0.5 * NA * probability_function(2);

%% Tx_1 transmit bit 0
mu_1 = mu_0 + NA * probability_function(1);

threshold = ceil((mu_1+mu_0)/log(mu_1/mu_0));
Pe = 0.5 * ( (1 - poisscdf(threshold, mu_0)) + poisscdf(threshold, mu_1));
end