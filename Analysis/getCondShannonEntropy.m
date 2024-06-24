function [H] = getCondShannonEntropy(xhit,xmiss,nBins,minVal,maxVal)

binEdges = linspace(minVal,maxVal,nBins+1);

% Binning the array 
binnedXHit = discretize(xhit,binEdges);
binnedXMiss = discretize(xmiss,binEdges);
binnedXAll = horzcat(binnedXHit,binnedXMiss);

nTrials = size(binnedXAll,2); % number of all trials
nHit = size(binnedXHit,2); % number of all hit trials 
nMiss = size(binnedXMiss,2); % number of all miss trials 

% Hits 
pkHit = zeros(nBins,1);
HkHit = zeros(nBins,1);
pHit = nHit/nTrials;
for k=1:nBins
    nkHit = sum(binnedXHit==k);
    pkHit(k) = nkHit/nHit;
    if pkHit(k) == 0
        HkHit(k) = 0;
    else
        HkHit(k) = -1*pkHit(k)*(log2(pkHit(k)));
    end
end
HXYHit = pHit*(sum(HkHit));
% Misses 
pkMiss = zeros(nBins,1);
HkMiss = zeros(nBins,1);
pMiss = nMiss/nTrials;
for k=1:nBins
    nkMiss = sum(binnedXMiss==k);
    pkMiss(k) = nkMiss/nMiss;
    if pkMiss(k) == 0
        HkMiss(k) = 0;
    else
        HkMiss(k) = -1*pkMiss(k)*(log2(pkMiss(k)));
    end
end
HXYMiss = pMiss*(sum(HkMiss));

H = HXYHit + HXYMiss;