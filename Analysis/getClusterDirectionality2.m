function [mode, fittedVmFinal] = getClusterDirectionality2(dirComb,r2thres,modeMax,nBins)

% % Description 

% Input: 
% Waves - wave structure
% dirComb - array of the directions 
% r2thres - threshold for r2 - godness of fit 
% modeMax - max number of peaks in the mix von misses fitting
% nBins - Number of bins for histogram fitting
% Output
% mode: 
% 0 -> non-directional 
% 1 -> unimodal
% 2 -> bimodal
% 3 -> multimodal
% fittedVm - the final fitted Vm 


% dirComb = horzcat(Waves(1:end).waveDir);
dirPDF.edges = linspace(-pi,pi,nBins+1);
[dirPDF.prob,dirPDF.edges,dirPDF.bins] = histcounts(dirComb,dirPDF.edges,'Normalization','probability');

% r2 = 0;
% k = 0;
% d = dirComb';
% d(d<0) = d(d<0)+2*pi;

r2 = zeros(1,modeMax);
for k = 1:modeMax
    fittedVm(k) = fitmvmdist(dirComb',k,'MaxIter',250,'ErrorThreshold',1e-5);
    fittedProb = zeros(1,size(dirPDF.prob,2));
    for i=1:size(dirPDF.prob,2)
        edgeStart = dirPDF.edges(i);
        edgeEnd = dirPDF.edges(i+1);
        phi = linspace(edgeStart,edgeEnd,20)';
        dphi = phi(2)-phi(1);
        fittedProb(i) = sum(fittedVm(k).pdf(phi),'all')*dphi;
    end
    r2(k) = r2score(sort(dirPDF.prob)',sort(fittedProb)');
    if r2(k)>r2thres
        r2max = r2(k);
        fittedVmFinal = fittedVm(k);
        break;
    end
    if k==modeMax
        r2max = r2(k);
        fittedVmFinal = fittedVm(k);
    end
end


% [p,~] =circ_rtest(dirComb);
% disp('Rayleigh test for non-uniformity');
% disp(['p-val = ',num2str(p)]);
if (k == 1 && any(fittedVmFinal.componentProportion < 0.2))
    mode = 0;
end
if (k == 1 && all(fittedVmFinal.componentProportion >= 0.2))
    mode = 1;
end
if (k == 2 && any(fittedVmFinal.componentProportion < 0.2))
    mode = 1;
end
if (k == 2 && all(fittedVmFinal.componentProportion >= 0.2))
    mode = 2;
end
if (k == 3 && any(fittedVmFinal.componentProportion < 0.2))
    mode = 2;
end
if (k == 3 && all(fittedVmFinal.componentProportion >= 0.2))
    mode = 3;
end

if mode ~= 0
    fittedVmFinal = fittedVm(mode);
end
disp(['Direction Mode: ' ,num2str(mode)]);

% 
% if p< 0.05 
%     
%     if k == 1
%         mode = 1;
%     end
%     if (k == 2 && any(fittedVmFinal.kappa < 2))
%         mode = 1;
%     end
%     if (k == 2 && all(fittedVmFinal.kappa >= 2))
%         mode = 2;
%     end
%     if (k == 3)
%         mode = 3;
%     end
% else
