function [probability_function] = probability_mobile(Rx1_estimate_coefficient, Rx2_estimate_coefficient, Rx3_estimate_coefficient, r, diffusion_coefficient, t, tau,distance)
% Description:
% r: radius of receiver
% distance: distance between transmitter and receiver
% tau : sample time in one time duration
% molecules_coefficient: molecule diffusion
% transmitter_coefficient: transmitter diffusion
% receiver_coefficient: receiver diffusion

initial_distance = distance;
V = 4*pi*r^3/3;
molecule_coefficient = diffusion_coefficient(1);
transmitter_coefficient = diffusion_coefficient(2);
receiver_coefficient = diffusion_coefficient(3);
D1 = molecule_coefficient + receiver_coefficient;
D2 = transmitter_coefficient + receiver_coefficient;
b = [Rx1_estimate_coefficient'; Rx2_estimate_coefficient'];
% h_1 : the probability of each molecules released by Tx1 observed by Rx1
% h_2 : the probability of each molecules released by Tx2 observed by Rx2
% h_3 : the probability of each molecules released by Tx3 observed by Rx3
h_1 = b(1,1)*V*erfc(initial_distance(1)/((4*D2)^b(1,2)*(D1/D2*tau+t)^b(1,3)))/(4*pi*D2*initial_distance(1));
h_2 = b(2,1)*V*erfc(initial_distance(2)/((4*D2)^b(2,2)*(D1/D2*tau+t)^b(2,3)))/(4*pi*D2*initial_distance(2));
h_3 = b(3,1)*V*erfc(initial_distance(3)/((4*D2)^b(3,2)*(D1/D2*tau+t)^b(3,3)))/(4*pi*D2*initial_distance(3));
probability_function = [h_1;h_2;h_3];
end