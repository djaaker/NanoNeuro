function [PGD] = getPGD(xgp,behaviourTrace,parameters)

if strcmp(parameters.experiment,'cue')
    totaltime = parameters.windowAfterCue + parameters.windowBeforeCue + parameters.ts;
elseif strcmp(parameters.experiment,'self')
    totaltime = parameters.windowAfterPull + parameters.windowBeforePull + parameters.ts;
end

time = [parameters.ts:parameters.ts:totaltime]; 

PGD = zeros(size(behaviourTrace,2),size(time,2));

for ii=1:size(behaviourTrace,2)
    p = xgp{1,ii};
    [pm,~,dx,dy] = phase_gradient_complex_multiplication(p,parameters.xspacing,parameters.yspacing );
    % Phase gradient directionality 
    PGD(ii,:) = phase_gradient_directionality(pm,dx,dy); 
end
  
end

