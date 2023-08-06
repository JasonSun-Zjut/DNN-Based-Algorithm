function [solution, freq_variable] = Diversification_Generation(variable_number, low_bound, upper_bound, sub_range_length, solution_number)
sub_range_number = (upper_bound - low_bound) / sub_range_length;
r = randi(sub_range_number, solution_number, variable_number);
[C, ~, ic] = unique(r);
a_counts = accumarray(ic, 1);
freq_variable = [C, a_counts];
solution = zeros(solution_number, variable_number);
for i = 1 : size(r, 1)
solution(i, 1) = randi([(r(i, 1) - 1) * sub_range_length + low_bound, r(i, 1) * sub_range_length + low_bound], 1);
solution(i, 2) = randi([(r(i, 2) - 1) * sub_range_length + low_bound, r(i, 2) * sub_range_length + low_bound], 1);
end
end 