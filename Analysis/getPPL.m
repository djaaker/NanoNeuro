function PPL = getPPL(IntanBehaviour,z_score,nIterrate,plotFlag,parameters,shank)

% Reference : Spontaneous travelling cortical waves gate perception in 
% behaving primates, Nature, 2020 
% shank = 1 - probe lfp is used 
% shank = 0 - gird lfp is used

if shank == 0
    nElectrodes = parameters.rows*parameters.cols;
else
    nElectrodes = parameters.nShank;
end

if shank == 0    
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    [PPL.PPLHit] = calPPL(xgp,parameters);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    [PPL.PPLMiss] =  calPPL(xgp,parameters);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.hitTrace, 'UniformOutput', false);
    [PPL.PPLHitReward] =  calPPL(xgp,parameters);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.missTrace, 'UniformOutput', false);
    [PPL.PPLFA] =  calPPL(xgp,parameters);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    [PPL.PPLMIHit] =  calPPL(xgp,parameters);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    [PPL.PPLMIFA] =  calPPL(xgp,parameters);
else
    xgp = arrayfun(@(s) s.xgpProbe, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    [PPL.PPLHit] =  calPPL(xgp,parameters);
    xgp = arrayfun(@(s) s.xgpProbe, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    [PPL.PPLMiss] =  calPPL(xgp,parameters);
    xgp = arrayfun(@(s) s.xgpProbe, IntanBehaviour.hitTrace, 'UniformOutput', false);
    [PPL.PPLHitReward] =  calPPL(xgp,parameters);
    xgp = arrayfun(@(s) s.xgpProbe, IntanBehaviour.missTrace, 'UniformOutput', false);
    [PPL.PPLFA] =  calPPL(xgp,parameters);
    xgp = arrayfun(@(s) s.xgpProbe, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    [PPL.PPLMIHit] =  calPPL(xgp,parameters);
    xgp = arrayfun(@(s) s.xgpProbe, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    [PPL.PPLMIFA] =  calPPL(xgp,parameters);
end

if plotFlag == 1
    figure();
    subplot(2,1,1);
    title("Percentage Phase across all electrodes - Hits")
    data = removeNaNRows(peakSort2DArray(reshape(PPL.PPLHit,[],size(PPL.PPLHit,3)),'descend',2));
    imagesc(IntanBehaviour.cueHitTrace(1).time,1:size(data,1),data); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','Cue','LabelVerticalAlignment','top');
    xline(mean(IntanBehaviour.reactionTime,'all'),'--w','Avg. Reaction Time','LabelVerticalAlignment','top');
    subplot(2,1,2);
    title("Percentage Phase across all electrodes - Misses")
    data = removeNaNRows(peakSort2DArray(reshape(PPL.PPLMiss,[],size(PPL.PPLMiss,3)),'descend',2));
    imagesc(IntanBehaviour.cueMissTrace(1).time,1:size(data,1),data); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','Cue','LabelVerticalAlignment','top');

    figure();
    subplot(2,1,1);
    title("Percentage Phase across all electrodes - Hits")
    imagesc(IntanBehaviour.hitTrace(1).time,1:nElectrodes,reshape(PPL.PPLHitReward,[],size(PPL.PPLHitReward,3))); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','Threshold','LabelVerticalAlignment','top');
    % yyaxis right; box off;
    % plot(IntanBehaviour.hitTrace(1).time,squeeze(mean(PPL.PPLHitReward,[1 2],'omitnan')),'-r','Linewidth',0.8);
    % xline(-mean(IntanBehaviour.reactionTime,'all'),'--w','Avg. Cue Time','LabelVerticalAlignment','top');
    subplot(2,1,2);
    title("Percentage Phase across all electrodes - FA")
    imagesc(IntanBehaviour.missTrace(1).time,1:nElectrodes,reshape(PPL.PPLFA,[],size(PPL.PPLFA,3))); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','Threshold','LabelVerticalAlignment','top');
    % yyaxis right; box off;
    % plot(IntanBehaviour.missTrace(1).time,squeeze(mean(PPL.PPLFA,[1 2],'omitnan')),'-r','Linewidth',0.8);

    figure();
    subplot(2,1,1);
    title("Percentage Phase across all electrodes - Hits MI")
    imagesc(IntanBehaviour.MIHitTrace(1).time,1:nElectrodes,reshape(PPL.PPLMIHit,[],size(PPL.PPLMIHit,3))); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','MI','LabelVerticalAlignment','top');
    subplot(2,1,2);
    title("Percentage Phase across all electrodes - Misses MI")
    imagesc(IntanBehaviour.MIFATrace(1).time,1:nElectrodes,reshape(PPL.PPLMIFA,[],size(PPL.PPLMIFA,3))); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','MI','LabelVerticalAlignment','top');


    figure();
    plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPL.PPLHit,[1 2])),'-r','LineWidth',1.2); hold on;
    plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPL.PPLMiss,[1 2])),'-k','LineWidth',0.2);
    ylabel("Percentage Phase Locking"); xlabel("Time (s)");
    xline(0,'--k','Cue','LabelVerticalAlignment','top');
    xline(mean(IntanBehaviour.reactionTime,'all'),'--k','Avg. Reaction Time','LabelVerticalAlignment','top');
    title('Percentage Phase Locking for hits');box off;legend('Hits','Misses');%,'Misses');

    figure();
    plot(IntanBehaviour.MIHitTrace(1).time,squeeze(nanmean(PPL.PPLMIHit,[1 2])),'-r','LineWidth',1.2); hold on;
    plot(IntanBehaviour.MIFATrace(1).time,squeeze(nanmean(PPL.PPLMIFA,[1 2])),'-k','LineWidth',0.2);
    ylabel("Percentage Phase Locking"); xlabel("Time (s)");
    xline(0,'--k','MI','LabelVerticalAlignment','top');
    title('Percentage Phase Locking for hits');box off;legend('Hits','FAs');%,'Misses');
end

if z_score == 1
    % z-scoring 
    xgpHit = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    xgpMiss = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    xgpComb = [xgpHit xgpMiss];
    nHit = size(xgpHit,2);
    nMiss = size(xgpMiss,2);
    nTot = nHit + nMiss;
    nullDistHit = zeros(parameters.rows,parameters.cols,size(PPL.PPLHit,3),nIterrate);
    nullDistMiss = zeros(parameters.rows,parameters.cols,size(PPL.PPLMiss,3),nIterrate);
    for j=1:nIterrate
        randIndex = randperm(nTot);
        xgpHitRand = xgpComb(randIndex(1:nHit));
        xgpMissRand = xgpComb(randIndex(nHit+1:end));
        nullDistHit(:,:,:,j) = calPPL(xgpHitRand,parameters);
        nullDistMiss(:,:,:,j) = calPPL(xgpMissRand,parameters);
        j
    end
    muHit = mean(nullDistHit,4); % Mean of the null distribution
    sigmaHit = std(nullDistHit,0,4); % Standard deviation of null distribution
    muMiss = mean(nullDistMiss,4); % Mean of the null distribution
    sigmaMiss = std(nullDistMiss,0,4); % Standard deviation of null distribution
    
    PPL.PPLHitz = (PPL.PPLHit-muHit)./sigmaHit;
    PPL.PPLMissz = (PPL.PPLMiss-muMiss)./sigmaMiss;
    
    figure();
    plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPL.PPLHitz,[1 2])),'-r','LineWidth',1.2); hold on;
    % plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAMiss,[1 2])),'-k','LineWidth',1);
    plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPL.PPLMissz,[1 2])),'-k','LineWidth',1);
    ylabel("z-score"); xlabel("Time (s)");
    xline(0,'--r','Cue','LabelVerticalAlignment','top');
    xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
    title('z-scored Phase Alignment for Hits');box off; legend('Hits','Miss');
    
    figure();
    subplot(2,1,1);
    title("Percentage Phase across all electrodes - Hits")
    imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PPL.PPLHitz,[],size(PPL.PPLHitz,3))); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','Cue','LabelVerticalAlignment','top');
    subplot(2,1,2);
    title("Percentage Phase across all electrodes - Misses")
    imagesc(IntanBehaviour.cueMissTrace(1).time,1:32,reshape(PPL.PPLMissz,[],size(PPL.PPLMissz,3))); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','Cue','LabelVerticalAlignment','top');
%     
%     nIterrate = 200;
%     xgpHit = arrayfun(@(s) s.xgp, IntanBehaviour.hitTrace, 'UniformOutput', false);
%     xgpMiss = arrayfun(@(s) s.xgp, IntanBehaviour.missTrace, 'UniformOutput', false);
%     xgpComb = [xgpHit xgpMiss];
%     nHit = size(xgpHit,2);
%     nMiss = size(xgpMiss,2);
%     nTot = nHit + nMiss;
%     nullDistHit = zeros(parameters.rows,parameters.cols,size(PPL.PPLHitReward,3),nIterrate);
%     nullDistMiss = zeros(parameters.rows,parameters.cols,size(PPL.PPLFA,3),nIterrate);
%     for j=1:nIterrate
%         randIndex = randperm(nTot);
%         xgpHitRand = xgpComb(randIndex(1:nHit));
%         xgpMissRand = xgpComb(randIndex(nHit+1:end));
%         nullDistHit(:,:,:,j) = calPPL(xgpHitRand,parameters);
%         nullDistMiss(:,:,:,j) = calPPL(xgpMissRand,parameters);
%         j
%     end
%     muHitReward = mean(nullDistHit,4); % Mean of the null distribution
%     sigmaHitReward = std(nullDistHit,0,4); % Standard deviation of null distribution
%     muFA = mean(nullDistMiss,4); % Mean of the null distribution
%     sigmaFA = std(nullDistMiss,0,4); % Standard deviation of null distribution
%     
%     PPL.PPLHitRewardz = (PPL.PPLHit-muHitReward)./sigmaHitReward;
%     PPL.PPLFAz = (PPL.PPLMiss-muFA)./sigmaFA;
end