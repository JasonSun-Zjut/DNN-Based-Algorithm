function Pe = error_probability(estimate_coefficient, n, NA, distance, diffusion_coefficient)
% Description: error probability of Rx in current time slot
% n: nth time slot
% NA: number of molecules to emit
r = 3;
probability_function = probability_mobile(estimate_coefficient(1), estimate_coefficient(2), estimate_coefficient(3), r, diffusion_coefficient, (n-1)*T, tau, distance);
mu_0 = 0.5 * NA * probability_function(2);

%% Tx_1 transmit bit 0
mu_1 = mu_0 + NA * probability_function(1);

threshold = ceil((mu_1+mu_0)/log(mu_1/mu_0));
Pe = 0.5 * ( (1 - poisscdf(threshold, mu_0)) + poisscdf(threshold, mu_1));
end