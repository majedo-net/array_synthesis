 function [element_pattern] = load_element(inputArg1)
%LOAD_ELEMENT Load element pattern CSV data from HFSS simulation
%   inputArg1 contains a string with filename
arr_data = table2array(readtable(inputArg1));
theta = linspace(0,180,181);
phi = transpose(linspace(0,360,361));
element_pattern = zeros(numel(phi),numel(theta));
for theta_d = 1:91
    idx1 = (theta_d-1)*180+theta_d;
    idx2 = idx1 + 180;
    element_pattern(1:181,theta_d) = interp1(0:2:181,arr_data(idx1+90:idx2,3),0:180);
    element_pattern(181:361,theta_d) = interp1(0:2:181,flip(arr_data(idx1:idx1+90,3)),0:180);
end

for phi_d = 1:361
    element_pattern(phi_d,:) = interp1(0:2:180,element_pattern(phi_d,1:91),0:180);
end
