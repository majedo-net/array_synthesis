function [samples] = rejectionSample(trueFunc,xs,Ns)
%REJ_SAMPLE Summary of this function goes here
%   Detailed explanation goes here
max_x = max(xs);
min_x = min(xs);
range_x = max_x - min_x;
max_f = max(trueFunc);
samples = [];
while numel(samples)<Ns
    u = rand(1).*range_x + min_x;
    v = rand(1).*max_f;
    fu = interp1(xs,trueFunc,u);
    if v <= fu
        samples(end+1) = u;
    end
end
end

