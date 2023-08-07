function best_solution = select_diversification_generation(input, trial_solution, estimate_coefficient)
% Description: Select the best solution among all neighborhood solutions
% Parameters
% input : the input of DNN
% trial_solution: neighborhood solutions quantized from the output of DNN
transmitter_diffusion_coefficient = zeros(3,1); % DTx = 10 × 10−13m2/s
receiver_diffusion_coefficient = zeros(1,1); % DRx = 10 × 10−13 m2 /s;
molecule_diffusion_coefficient = zeros(1,1); % DA = 5 × 10−9m2/s
n = 5;
diffusion_coefficient = [transmitter_diffusion_coefficient; receiver_diffusion_coefficient; molecule_diffusion_coefficient];
best_solution = zeros(1,3);
min_Pe = realmax;
% select the best solution 
for i = 1 : 1 : size(trial_solution)
    % calculate the BER of each trial solution
    Pe = error_probability(estimate_coefficient, n, trial_solution(i,:), input, diffusion_coefficient);
    if(Pe < min_Pe)
        min_Pe = Pe;
    end
    best_solution = trial_solution(i, :);
end