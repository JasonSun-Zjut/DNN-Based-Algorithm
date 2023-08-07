function [coordinate] = rand_coordinate_generate(point, distance)

% Description:
% find a point that is a certain distance from another point
% Parameter:

% x-axis y-axis 扩展的步长
%{
step = distance;
x = point(1)  : point(1) + step;
y = point(2)  : point(2) + step;
[X,Y] = meshgrid(x,y);
distance_vector = sqrt((X-point(1)).^2 + (Y-point(2)).^2);
[row, col] = find(distance_vector == distance);
% 符合条件的坐标
coordiante_vector = [x(row)',y(col)'];
%}
% find a point that is a certain distance from another point
% x_coordinate is equal to y_coordinate
x_coordinate = sqrt(distance^2/2);
y_coordinate = sqrt(distance^2/2);
coordinate = [point(1) + x_coordinate, point(2) + y_coordinate];

end