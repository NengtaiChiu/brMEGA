function [trend]=median_filter(X,K)
%X = raw signal
%==========================================================================
%prepare data
trend = zeros(size(X));
nX = [zeros(1000,1) ;X ;zeros(1000,1)];

%%
for i = 1001:(length(X)+1000)
    m = median(nX(i-K:i+K));
    trend(i-1000) = m;
end


 





