function [d] = RaisedPowerSeries(fmax,r,N)
%%
% @file RaisedPowerSeries.m
%
% @brief Generate raised power series vector for a given max frequency
% 
% input:
%   fmax: Frequency in hertz to evaluate array factor
%   
%   r: raised power exponent
%
%   N: Quadrant number of elements. Total size for a 2d planar array will
%   be (2N+1)^2 
% 
%
% @copyright Copyright (c) 2022 Augustus Aerospace Company, all rights reserved.
%
if nargin < 3
    disp('Missing arguments');
    return;
end
d = zeros(N,1);
d0 = (3e8/fmax)/2;
xi = 1/((N^r)-(N-1)^r);
for n=1:N
    d(n) = d0*xi*n^r;
end
d = cat(1,d(1),d);
end