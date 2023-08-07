% Description:
% The DNN-Based Optimization Algorithm
clear;
% the number of transmitters in this model
node = 3;
% init a neural network with three hidden layers
net = feedforwardnet([2 * node + 2, 2 * node, 2 * node]);
net.trainfcn='trainbr';
net.layers{2}.transferFcn='softmax';
net.layers{3}.transferFcn='poslin';

dataset_number = 1000;

batch_train_dataset_size = 100;
% dataset generation
[dataset, estimate_coefficient] = dataset_generation(node, dataset_number);

% train network
net = train(net, dataset(1:dataset_number, 1:node)', dataset(1:dataset_number, node + 1:end)');

for i = 1 : size(dataset, 1)
     net_output = sim(net, dataset(i,1 : node)');
     sub_range_length = 100;
     K = 10;
     [trial_solution, ~] = Diversification_Generation(node, molecule_allocate_low_bound, molecule_allocate_upper_bound, sub_range_length, quantization_data_number);
     best_solution = select_diversification_generation(dataset(i, 1:node), trial_solution, estimate_coefficient);
     % queue init
     queue_memory = Queue;
     queue_memory.enqueue([dataset(i,1 : node), best_solution]);
     if mod (i, batch_train_dataset_size) == 0
         retrain_data_set = queue_memory.dequeue;
     end
     % retrain
     net = train(net, retrain_data_set(1 : batch_train_dataset_size, 1:node)', dataset(1 : batch_train_dataset_size, node + 1:end)');
end