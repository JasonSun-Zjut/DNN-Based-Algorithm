function [probability_function] = probability_mobile(Rx1_estimate_coefficient, Rx2_estimate_coefficient, r, diffusion_coefficient, t, tau,distance)
% Description:
% r: radius of receiver
% distance: distance between transmitter and receiver
% tau : sample time in one time duration
% molecules_coefficient: molecule diffusion
% transmitter_coefficient: transmitter diffusion
% receiver_coefficient: receiver diffusion

% only the initial distance between transmitter and receiver is known

% mobile fitted coefficient [500 10 10] topology [d = 10 r = 4 h = 2]
initial_distance = distance;
V = 4*pi*r^3/3;
molecule_coefficient = diffusion_coefficient(1);
transmitter_coefficient = diffusion_coefficient(2);
receiver_coefficient = diffusion_coefficient(3);
D1 = molecule_coefficient + receiver_coefficient;
D2 = transmitter_coefficient + receiver_coefficient;
b = [Rx1_estimate_coefficient'; Rx2_estimate_coefficient'];
h_1 = b(1,1)*V*erfc(initial_distance(1)/((4*D2)^b(1,2)*(D1/D2*tau+t)^b(1,3)))/(4*pi*D2*initial_distance(1));
h_2 = b(2,1)*V*erfc(initial_distance(2)/((4*D2)^b(2,2)*(D1/D2*tau+t)^b(2,3)))/(4*pi*D2*initial_distance(2));
%h = V*exp(-initial_distance^2/(4*(D1*tau+D2*t)))/(4*pi*(D1*tau+D2*t))^(3/2);
probability_function = [h_1;h_2];
end