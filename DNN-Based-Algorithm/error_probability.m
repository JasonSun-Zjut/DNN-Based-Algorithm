function Pe = error_probability(estimate_coefficient, n, NA, distance, diffusion_coefficient)
% Description: error probability of Rx in current time slot
% n: nth time slot
% NA: number of molecules to emit

Pe = zeros(3,1);
mu_0 = zeros(3,1);
mu_1 = zeros(3,1);
r = 3;
probability_function = probability_mobile(estimate_coefficient(1), estimate_coefficient(2), estimate_coefficient(3), r, diffusion_coefficient, (n-1)*T, tau, distance);
% Tx_1 transmit bit 0
mu_0(1) = 0.5 * NA * probability_function(1);
% Tx_1 transmit bit 1
mu_1(1) = mu_0(1) + NA * probability_function(1);
threshold = ceil((mu_1(1) + mu_0(1))/log(mu_1(1) / mu_0(1)));
Pe(1) = 0.5 * ( (1 - poisscdf(threshold, mu_0(1))) + poisscdf(threshold, mu_1(1)));

% Tx_2 transmit bit 0
mu_0(2) = 0.5 * NA * probability_function(2);
% Tx_2 transmit bit 1
mu_1(2) = mu_0(2) + NA * probability_function(2);
threshold = ceil((mu_1(2) + mu_0(2))/log(mu_1(2) / mu_0(2)));
Pe(2) = 0.5 * ( (1 - poisscdf(threshold, mu_0(2))) + poisscdf(threshold, mu_1(2)));

% Tx_3 transmit bit 0
mu_0(3) = 0.5 * NA * probability_function(3);
% Tx_3 transmit bit 1
mu_1(3) = mu_0(3) + NA * probability_function(3);
threshold = ceil((mu_1(3) + mu_0(3))/log(mu_1(3) / mu_0(3)));
Pe(3) = 0.5 * ( (1 - poisscdf(threshold, mu_0(3))) + poisscdf(threshold, mu_1(3)));

end