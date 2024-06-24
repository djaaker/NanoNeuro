function [MI] = getMutualInformation(xHit,xMiss,parameters)
% Reference : Propagating waves mediate information transfer in the motor
% cortex, Nature Neuro, 2006 

nHits = size(xHit,2);
nMiss = size(xMiss,2);

Nbins = ceil(exp(0.626 + (0.4*log(nHits+nMiss-1)))); % number of bins
% Ref - Comparison of Hilbert transform and wavelet methods for the
%       analysis of neuronal synchrony, JNeuroMethods, 2001

HX = zeros(parameters.rows,parameters.cols,size(xHit{1,1},3)); % Shannon entropy 
HXY = zeros(parameters.rows,parameters.cols,size(xHit{1,1},3)); % conditonal shannon entropy 
MI = zeros(parameters.rows,parameters.cols,size(xHit{1,1},3)); % Mutual information 

for t=1:size(xHit{1,1},3)
    for i=1:parameters.rows
        for j=1:parameters.cols
            xijtHit = cellfun(@(s) s(i,j,t),xHit);
            xijtMiss = cellfun(@(s) s(i,j,t),xMiss);
            xijtAll = horzcat(xijtHit,xijtMiss);
            if sum(isnan(xijtAll))>0
                HX(i,j,t) = NaN;
                HXY(i,j,t) = NaN;
                MI(i,j,t) = NaN;
            else
                HX(i,j,t) = getShannonEntropy(xijtAll,Nbins,-pi,pi);
                HXY(i,j,t) = getCondShannonEntropy(xijtHit,xijtMiss,Nbins,-pi,pi);
                MI(i,j,t) = HX(i,j,t) - HXY(i,j,t);
            end
        end
    end
end
