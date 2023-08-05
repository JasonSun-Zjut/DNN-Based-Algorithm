%%
clear;
node = 2;
net = feedforwardnet([2 * node + 2, 2 * node, 2 * node]);
net.trainfcn='trainbr';
net.layers{2}.transferFcn='softmax';
net.layers{3}.transferFcn='poslin';
molecule_allocate_low_bound = 5000;
molecule_allocate_upper_bound = 10000;
distance_low_bound = 10;
distance_upper_bound = 30;
Init_Data_Number = 20;
Rand_Input = sort(randi([distance_low_bound, distance_upper_bound], Init_Data_Number, node), 2);
Rand_Output = randi([molecule_allocate_low_bound, molecule_allocate_upper_bound], 20, node);
Rand_Data_Set = [Rand_Input, Rand_Output];
% init network
net = train(net, Rand_Input', Rand_Output');

% train network
Train_Data_Number = 100;
Train_Input = sort(randi([distance_low_bound, distance_upper_bound], Train_Data_Number, node), 2);
for i = 1 : size(Train_Input, 1)
     net_output = sim(net, Train_Input(i,:)');
     sub_range_length = 100;
     [trail_solution, freq_variable] = Diversification_Generation(node, molecule_allocate_low_bound, molecule_allocate_upper_bound, sub_range_length, Train_Data_Number);
     % error probabiltiy calculate
end