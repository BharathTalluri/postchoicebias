function y=nanlength(x)

%function y=nanlength(x)
%returns length of x without nans

y=length(find(~isnan(x)));