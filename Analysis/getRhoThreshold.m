function [rhoThres] = getRhoThreshold(behaviourTrace1,behaviourTrace2,parameters,nShuffle,threshold,fraction)

% warningid = 'MATLAB:smoothn:SUpperBound';
warning('off','all');

X = parameters.X;
Y = parameters.Y;

nTrialShuffle1 = ceil(fraction*size(behaviourTrace1,2));
a1 = randperm(size(behaviourTrace1,2));
trials1 = a1(1:nTrialShuffle1);
nTrialShuffle2 = ceil(fraction*size(behaviourTrace2,2));
a2 = randperm(size(behaviourTrace2,2));
trials2 = a2(1:nTrialShuffle2);

rhoAllShuffle = [];
rhoAll = [];

for i=1:nTrialShuffle1
    evaluationPoints = find_evaluation_points(behaviourTrace1(trials1(i)).xgp,0,0.2);
    for j=1:nShuffle
        if j==1
            pShuffle = behaviourTrace1(trials1(i)).xgp;
%             for ll=1:size(pShuffle,3)
%                 pShuffle1(:,:,ll) = inpaint_nans(pShuffle(:,:,ll),3);
%             end
            [~,~,dx,dy] = phase_gradient_complex_multiplication( pShuffle, parameters.xspacing,parameters.yspacing );
            source = find_source_points( evaluationPoints, X, Y, dx, dy );
            for jj = 1:length(evaluationPoints)
                ph = angle( pShuffle(:,:,evaluationPoints(jj)) );
                rho(jj) = phase_correlation_distance( ph, source(:,jj), parameters.xspacing,parameters.yspacing );
            end
            rhoAll = [rhoAll,rho];
        else
            pShuffle = shuffle_channels(behaviourTrace1(trials1(i)).xgp);
%             for ll=1:size(pShuffle,3)
%                 pShuffle1(:,:,ll) = inpaint_nans(pShuffle(:,:,ll),3);
%             end
            [~,~,dx,dy] = phase_gradient_complex_multiplication( pShuffle, parameters.xspacing,parameters.yspacing );
            source = find_source_points( evaluationPoints, X, Y, dx, dy );
            for jj = 1:length(evaluationPoints)
                ph = angle( pShuffle(:,:,evaluationPoints(jj)) );
                rho(jj) = phase_correlation_distance( ph, source(:,jj), parameters.xspacing,parameters.yspacing );
            end
            rhoAllShuffle = [rhoAllShuffle,rho];
        end
    end
    disp([num2str(i) '/' num2str(nTrialShuffle1)]);
end

for i=1:nTrialShuffle2
    evaluationPoints = find_evaluation_points(behaviourTrace2(trials2(i)).xgp,0,0.2);
    for j=1:nShuffle
        if j==1
            pShuffle = behaviourTrace2(trials2(i)).xgp;
%             for ll=1:size(pShuffle,3)
%                 pShuffle1(:,:,ll) = inpaint_nans(pShuffle(:,:,ll),3);
%             end
            [~,~,dx,dy] = phase_gradient_complex_multiplication( pShuffle, parameters.xspacing,parameters.yspacing );
            source = find_source_points( evaluationPoints, X, Y, dx, dy );
            for jj = 1:length(evaluationPoints)
                ph = angle( pShuffle(:,:,evaluationPoints(jj)) );
                rho(jj) = phase_correlation_distance( ph, source(:,jj), parameters.xspacing,parameters.yspacing );
            end
            rhoAll = [rhoAll,rho];
        else
            pShuffle = shuffle_channels(behaviourTrace2(trials2(i)).xgp);
%             for ll=1:size(pShuffle,3)
%                 pShuffle1(:,:,ll) = inpaint_nans(pShuffle(:,:,ll),3);
%             end
            [~,~,dx,dy] = phase_gradient_complex_multiplication( pShuffle, parameters.xspacing,parameters.yspacing );
            source = find_source_points( evaluationPoints, X, Y, dx, dy );
            for jj = 1:length(evaluationPoints)
                ph = angle( pShuffle(:,:,evaluationPoints(jj)) );
                rho(jj) = phase_correlation_distance( ph, source(:,jj), parameters.xspacing,parameters.yspacing );
            end
            rhoAllShuffle = [rhoAllShuffle,rho];
        end
    end
    disp([num2str(i) '/' num2str(nTrialShuffle2)]);
end


% Plotting rho distribution

% Should the threshold calculated with all the rho (shuffled+unshuffled) or
% just with the shuffled rhos
rhoThres = prctile([rhoAllShuffle,rhoAll],threshold,'all');
% rhoThres = prctile([rhoAllShuffle],threshold,'all');
figure();histogram(rhoAll,'FaceColor','r'); hold on;
histogram(rhoAllShuffle,'FaceColor','b'); 
xline(rhoThres,'-r',{num2str(threshold) ' Percentile'});
xlabel('rho');
ylabel('Frequency');
title('Histogram of rho values for spatial shuffled electrodes');
drawnow
warning('on','all');
end
