function [indices_a, indices_b] = stratified_distribution(array_a, array_b, numbins)
% adapted from JW's python code
% Finds the common underlying distribution

if ~exist('numbins','var')
    numbins = 25; % default value
end

minimum = min(min(array_a),min(array_b));
maximum = max(max(array_a),max(array_b));
full_range = maximum - minimum;
binedges = linspace(minimum-(full_range/20), maximum+(full_range/20), numbins+1);
indices_a = 1:length(array_a);
indices_b = 1:length(array_b);
for i = 1:numbins
    % get the indices of values in both arrays that fall in this bin
    ind_a = find(array_a >= binedges(i) & array_a < binedges(i+1));
    ind_b = find(array_b >= binedges(i) & array_b < binedges(i+1));
    if length(ind_a) < length(ind_b)
        indices_b(datasample(ind_b,length(ind_b) - length(ind_a),'Replace',false)) = NaN;
    elseif length(ind_a) > length(ind_b)
        indices_a(datasample(ind_a,length(ind_a) - length(ind_b),'Replace',false)) = NaN;
    end
end
indices_a(isnan(indices_a)) = [];
indices_b(isnan(indices_b)) = [];