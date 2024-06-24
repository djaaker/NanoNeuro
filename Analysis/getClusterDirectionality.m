function [mode, fittedVm] = getClusterDirectionality(Waves,r2thres,modeMax)

% % Description 

% Input: 
% Waves - wave structure 
% r2thres - threshold for r2 - godness of fit 
% modeMax - max number of peaks in the mix von misses fitting
% 
% Output
% mode: 
% 0 -> non-directional 
% 1 -> unimodal
% 2 -> bimodal
% 3 -> multimodal
% fittedVm - the final fitted Vm 


dirComb = horzcat(Waves(1:end).waveDir);

r2 = 0;
k = 0;
d = dirComb';
d(d<0) = d(d<0)+2*pi;

while r2<r2thres
    k=k+1;
    fittedVm = fitmvmdist(dirComb',k,'MaxIter',250,'ErrorThreshold',1e-5);
%     fittedVm.kappa
    fittedData = fittedVm.random(size(d,1),false);
    fittedData(fittedData<0) = fittedData(fittedData<0)+2*pi;
    r2 = r2score(sort(d),sort(fittedData));
    if r2>r2thres
        disp(num2str(r2))
        break
    end
    if k == modeMax
        disp(['Max Mode condition reached with r2 score of ', num2str(r2)]);
        break
    end
end

[p,~] =circ_rtest(dirComb);
disp('Rayleigh test for non-uniformity');
disp(['p-val = ',num2str(p)]);

if p< 0.05 
    
    if k == 1
        mode = 1;
    end
    if (k == 2 && any(fittedVm.kappa < 2))
        mode = 1;
    end
    if (k == 2 && all(fittedVm.kappa >= 2))
        mode = 2;
    end
    if (k == 3)
        mode = 3;
    end
else
    if (k == 1 && any(fittedVm.kappa < 2))
        mode = 0;
    end
    if (k == 1 && all(fittedVm.kappa >= 2))
        mode = 1;
    end
    if (k == 2 && any(fittedVm.kappa < 2))
        mode = 1;
    end
    if (k == 2 && all(fittedVm.kappa >= 2))
        mode = 2;
    end
    if (k == 3 && any(fittedVm.kappa < 2))
        mode = 2;
    end
    if (k == 3 && all(fittedVm.kappa >= 2))
        mode = 3;
    end
end

disp(['Direction Mode: ' ,num2str(mode)]);
