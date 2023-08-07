function [dataset, estimate_coefficient] = dataset_generation(node, dataset_number)
% Description: 
% generate dataset based on the simulation result obtained by diffusion_rx_runner 
r = 3;
h = 1;
transmitter_diffusion_coefficient = zeros(3,1); % DTx = 10 × 10−13m2/s
receiver_diffusion_coefficient = zeros(1,1); % DRx = 10 × 10−13 m2 /s;
molecule_diffusion_coefficient = zeros(1,1); % DA = 5 × 10−9m2/s
initial_distance = zeros(3,1); % d011 = 15μm, d021 = 18μm, d031 = 25μm;
molecule_bound = zeros(2,1); %[molecule_low_boud, molecule_upper_boud]
distance_bound = zeros(2,1); %[distance_low_boud, distance_upper_boud]
dataset = zeros(dataset_number, node + node); %[d011, d021, d031, N1, N2, N3]
dataset_input = distance_bound(1) + (distance_bound(2) - distance_bound(1)) * rand(dataset_number, node);
dataset(1:dataset_number, 1:node) = dataset_input;

for i = 1 : 1 : dataset_number
tx_1_point = [0, 0, 0];
tx_2_point = [0, 0, -(2 * r + h)];
tx_3_point = [0, 0, -2 * (2 * r + h)];
rx_1_point = rand_coordinate_generate(tx_1_point, dataset_input(i, 1));
rx_2_point = rand_coordinate_generate(tx_2_point, dataset_input(i, 2));
rx_3_point = rand_coordinate_generate(tx_3_point, dataset_input(i, 3));

coordinate_vector = [tx_1_point; tx_2_point; tx_3_point; rx_1_point; rx_2_point; rx_3_point];
diffusion_coefficient = [transmitter_diffusion_coefficient; receiver_diffusion_coefficient; molecule_diffusion_coefficient];
NA = (molecule_bound(1) + molecule_bound(2)) / 2;
[res_Signal, b_Rx_1, b_Rx_2, b_Rx_3] = diffusion_rx_runner(NA, coordinate_vector, diffusion_coefficient);
dataset_output = zeros(1,3);
dataset_output(1,1) = sum(res_Signal(1, :));
dataset_output(1,2) = sum(res_Signal(2, :));
dataset_output(1,3) = sum(res_Signal(3, :));
dataset(i, node + 1 : end) = dataset_output;
estimate_coefficient = [b_Rx_1, b_Rx_2, b_Rx_3];
end
end

