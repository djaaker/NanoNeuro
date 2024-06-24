function [PA] = getPA(IntanBehaviour,z_score,nPerm,plotFlag,parameters,shank)

% Reference : Spontaneous travelling cortical waves gate perception in 
% behaving primates, Nature, 2020 
% shank = 1 - probe lfp is used 
% shank = 0 - gird lfp is used

if shank == 0
    nElectrodes = parameters.rows*parameters.cols;
else
    nElectrodes = parameters.nShank;
end
if strcmp(parameters.experiment,'cue')
    if shank == 0
        xgpHit = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
        xgpMiss = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    else
        xgpHit = arrayfun(@(s) s.xgpProbe, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
        xgpMiss = arrayfun(@(s) s.xgpProbe, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    end
    PA.Hit = calPhaseAlignment(xgpHit,parameters);
    PA.Miss = calPhaseAlignment(xgpMiss,parameters);
    [~,PA.PAPeakHit] = max(mean(PA.Hit,[1 2],'omitnan'));
    PA.HitAngles = rmmissing(reshape(cell2mat(cellfun(@(s) reshape(angle(s(:,:,PA.PAPeakHit)),nElectrodes,1), xgpHit,'UniformOutput',false)),[],1));
    [~,PA.PAPeakMiss] = max(mean(PA.Miss,[1 2],'omitnan'));
    PA.MissAngles = rmmissing(reshape(cell2mat(cellfun(@(s) reshape(angle(s(:,:,PA.PAPeakMiss)),nElectrodes,1), xgpMiss,'UniformOutput',false)),[],1));
end

if shank == 0
    xgpHitReward = arrayfun(@(s) s.xgp, IntanBehaviour.hitTrace, 'UniformOutput', false);
    xgpFA = arrayfun(@(s) s.xgp, IntanBehaviour.missTrace, 'UniformOutput', false);
else
    xgpHitReward = arrayfun(@(s) s.xgpProbe, IntanBehaviour.hitTrace, 'UniformOutput', false);
    xgpFA = arrayfun(@(s) s.xgpProbe, IntanBehaviour.missTrace, 'UniformOutput', false);
end
PA.HitReward = calPhaseAlignment(xgpHitReward,parameters);
PA.FA = calPhaseAlignment(xgpFA,parameters);
[~,PA.PAPeakHitReward] = max(mean(PA.HitReward,[1 2],'omitnan'));
PA.HitRewardAngles = rmmissing(reshape(cell2mat(cellfun(@(s) reshape(angle(s(:,:,PA.PAPeakHitReward)),nElectrodes,1), xgpHitReward,'UniformOutput',false)),[],1));
[~,PA.PAPeakFA] = max(mean(PA.FA,[1 2],'omitnan'));
PA.FAAngles = rmmissing(reshape(cell2mat(cellfun(@(s) reshape(angle(s(:,:,PA.PAPeakFA)),nElectrodes,1), xgpFA,'UniformOutput',false)),[],1));


if shank == 0
    xgpMIHit = arrayfun(@(s) s.xgp, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    xgpMIFA = arrayfun(@(s) s.xgp, IntanBehaviour.MIFATrace, 'UniformOutput', false);
else
    xgpMIHit = arrayfun(@(s) s.xgpProbe, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    xgpMIFA = arrayfun(@(s) s.xgpProbe, IntanBehaviour.MIFATrace, 'UniformOutput', false);
end
PA.MIHit = calPhaseAlignment(xgpMIHit,parameters);
PA.MIFA = calPhaseAlignment(xgpMIFA,parameters);
[~,PA.PAPeakMIHit] = max(mean(PA.MIHit,[1 2],'omitnan'));
PA.MIHitAngles = rmmissing(reshape(cell2mat(cellfun(@(s) reshape(angle(s(:,:,PA.PAPeakMIHit)),nElectrodes,1), xgpMIHit,'UniformOutput',false)),[],1));
[~,PA.PAPeakMIFA] = max(mean(PA.MIFA,[1 2],'omitnan'));
PA.MIFAAngles = rmmissing(reshape(cell2mat(cellfun(@(s) reshape(angle(s(:,:,PA.PAPeakMIFA)),nElectrodes,1), xgpMIFA,'UniformOutput',false)),[],1));


if plotFlag == 1
    if strcmp(parameters.experiment,'cue')
        figure();
        title("Phase Alignment across all electrodes - Hits")
        subplot(2,1,1);
        imagesc(IntanBehaviour.cueHitTrace(1).time,1:nElectrodes,peakSort2DArray(reshape(PA.Hit,[],size(PA.Hit,3)),'descend',2));  colormap(hot);
        ylabel("Electrodes");xlabel("Time (s)");
        xline(0,'-w','Cue','LabelVerticalAlignment','top');
        subplot(2,1,2);
        title("Phase Alignment across all electrodes - Misses")
        imagesc(IntanBehaviour.cueMissTrace(1).time,1:nElectrodes,peakSort2DArray(reshape(PA.Miss,[],size(PA.Miss,3)),'descend',2));  colormap(hot);
        ylabel("Electrodes");xlabel("Time (s)");
        xline(0,'-w','Cue','LabelVerticalAlignment','top');
    end
    
    figure();
    title("Phase Alignment across all electrodes - Hit Rewards")
    subplot(2,1,1);
    imagesc(IntanBehaviour.hitTrace(1).time,1:nElectrodes,reshape(PA.HitReward,[],size(PA.HitReward,3))); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','Reward','LabelVerticalAlignment','top');
    subplot(2,1,2);
    title("Phase Alignment across all electrodes - FA")
    imagesc(IntanBehaviour.missTrace(1).time,1:nElectrodes,reshape(PA.FA,[],size(PA.FA,3))); colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','No Reward','LabelVerticalAlignment','top');

    figure();
    title("Phase Alignment across all electrodes - Hit MI")
    subplot(2,1,1);
    imagesc(IntanBehaviour.MIHitTrace(1).time,1:nElectrodes,peakSort2DArray(reshape(PA.MIHit,[],size(PA.MIHit,3)),'descend',2));  colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','MI','LabelVerticalAlignment','top');
    subplot(2,1,2);
    title("Phase Alignment across all electrodes - FA MI")
    imagesc(IntanBehaviour.MIFATrace(1).time,1:nElectrodes,peakSort2DArray(reshape(PA.MIFA,[],size(PA.MIFA,3)),'descend',2));colormap(hot);
    ylabel("Electrodes");xlabel("Time (s)");
    xline(0,'-w','MI','LabelVerticalAlignment','top');
    
    if strcmp(parameters.experiment,'cue')
        figure();
        title("Phase Alignment averaged across Electrodes")
        plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PA.Hit,[1 2])),'-r','LineWidth',1.2); hold on;
        plot(IntanBehaviour.cueMissTrace(1).time,squeeze(nanmean(PA.Miss,[1 2])),'-k','LineWidth',1);
        % plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAFA,[1 2])),'-k','LineWidth',1);
        ylabel("Phase Alignment"); xlabel("Time (s)");
        xline(0,'--r','Cue','LabelVerticalAlignment','top');
        xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
        title('Phase Alignment for Hits');box off;legend('Hits','Miss');
    end

    figure();
    title("Phase Alignment averaged across Electrodes")
    plot(IntanBehaviour.hitTrace(1).time,squeeze(nanmean(PA.HitReward,[1 2])),'-r','LineWidth',1.2); hold on;
    plot(IntanBehaviour.missTrace(1).time,squeeze(nanmean(PA.FA,[1 2])),'-k','LineWidth',1);
    % plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAFA,[1 2])),'-k','LineWidth',1);
    ylabel("Phase Alignment"); xlabel("Time (s)");
    xline(0,'--r','Reward','LabelVerticalAlignment','top');
    xline(-mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Cue Time','LabelVerticalAlignment','top');
    title('Phase Alignment for Hits');box off;legend('Hit Reward','FA');

    figure();
    title("Phase Alignment averaged across Electrodes")
    plot(IntanBehaviour.MIHitTrace(1).time,squeeze(nanmean(PA.MIHit,[1 2])),'-r','LineWidth',1.2); hold on;
    plot(IntanBehaviour.MIFATrace(1).time,squeeze(nanmean(PA.MIFA,[1 2])),'-k','LineWidth',1);
    ylabel("Phase Alignment"); xlabel("Time (s)");
    xline(0,'--r','MI','LabelVerticalAlignment','top');
    title('Phase Alignment for Hits vs FA - Motion Initiation');box off;legend('Hits','FAs');
end

% z-scoring
if z_score == 1
    if strcmp(parameters.experiment,'cue')
        xgpComb = [xgpHit xgpMiss];
        nHit = size(xgpHit,2);
        nMiss = size(xgpMiss,2);
        nTot = nHit + nMiss;
        nullDistHit = zeros(size(PA.Hit,1),size(PA.Hit,2),size(PA.Hit,3),nPerm);
        nullDistMiss = zeros(size(PA.Miss,1),size(PA.Miss,2),size(PA.Miss,3),nPerm);
        for j=1:nPerm
            randIndex = randperm(nTot);
            xgpHitRand = xgpComb(randIndex(1:nHit));
            xgpMissRand = xgpComb(randIndex(nHit+1:end));
            nullDistHit(:,:,:,j) = calPhaseAlignment(xgpHitRand,parameters);
            nullDistMiss(:,:,:,j) = calPhaseAlignment(xgpMissRand,parameters);
            j
        end
        PA.muHit = mean(nullDistHit,4); % Mean of the null distribution
        PA.sigmaHit = std(nullDistHit,0,4); % Standard deviation of null distribution
        PA.muMiss = mean(nullDistMiss,4); % Mean of the null distribution
        PA.sigmaMiss = std(nullDistMiss,0,4); % Standard deviation of null distribution
        
        PA.Hitz = (PA.Hit-PA.muHit)./PA.sigmaHit;
        PA.Missz = (PA.Miss-PA.muMiss)./PA.sigmaMiss;
        if plotFlag == 1
            figure();
            plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PA.Hitz,[1 2])),'-r','LineWidth',1.2); hold on;
            plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PA.Missz,[1 2])),'-k','LineWidth',1);
            ylabel("z-score"); xlabel("Time (s)");
            xline(0,'--r','Cue','LabelVerticalAlignment','top');
            xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
            yline(2.56);        
            title('z-scored Phase Alignment averaged across electrodes');box off; legend('Hits','Miss');  
            
            figure();
            subplot(2,1,1);
            title("Phase Alignment across all electrodes - Hits")
            imagesc(IntanBehaviour.cueHitTrace(1).time,1:nElectrodes,reshape(PA.Hitz,[],size(PA.Hitz,3))); colormap(hot);
            ylabel("Electrodes");xlabel("Time (s)");
            xline(0,'-w','Cue','LabelVerticalAlignment','top');caxis([-3 3]);
            xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
            subplot(2,1,2);
            title("Phase Alignment across all electrodes - Misses")
            imagesc(IntanBehaviour.cueMissTrace(1).time,1:nElectrodes,reshape(PA.Missz,[],size(PA.Missz,3))); colormap(hot);
            ylabel("Electrodes");xlabel("Time (s)");
            xline(0,'-w','Cue','LabelVerticalAlignment','top');caxis([-3 3]);
        end
    end

    xgpComb = [xgpHitReward xgpFA];
    nHit = size(xgpHitReward,2);
    nMiss = size(xgpFA,2);
    nTot = nHit + nMiss;

    nullDistHit = zeros(size(PA.HitReward,1),size(PA.HitReward,2),size(PA.HitReward,3),nPerm);
    nullDistMiss = zeros(size(PA.FA,1),size(PA.FA,2),size(PA.FA,3),nPerm);
    for j=1:nPerm
        randIndex = randperm(nTot);
        xgpHitRand = xgpComb(randIndex(1:nHit));
        xgpMissRand = xgpComb(randIndex(nHit+1:end));
        nullDistHit(:,:,:,j) = calPhaseAlignment(xgpHitRand,parameters);
        nullDistMiss(:,:,:,j) = calPhaseAlignment(xgpMissRand,parameters);
        j
    end
    PA.muHitReward = mean(nullDistHit,4); % Mean of the null distribution
    PA.sigmaHitReward = std(nullDistHit,0,4); % Standard deviation of null distribution
    PA.muFA = mean(nullDistMiss,4); % Mean of the null distribution
    PA.sigmaFA = std(nullDistMiss,0,4); % Standard deviation of null distribution
    
    PA.HitRewardz = (PA.HitReward-PA.muHitReward)./PA.sigmaHitReward;
    PA.FAz = (PA.FA-PA.muFA)./PA.sigmaFA;
    if plotFlag == 1
        figure();
        plot(IntanBehaviour.hitTrace(1).time,squeeze(nanmean(PA.HitRewardz,[1 2])),'-r','LineWidth',1.2); hold on;
        plot(IntanBehaviour.missTrace(1).time,squeeze(nanmean(PA.FAz,[1 2])),'-k','LineWidth',1);
        ylabel("z-score"); xlabel("Time (s)");
        xline(0,'--r','Reward','LabelVerticalAlignment','top');
        xline(-mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Cue Time','LabelVerticalAlignment','top');
        title('z-scored Phase Alignment averaged across electrodes');box off; legend('Hits','FA');  
        
        figure();
        subplot(2,1,1);
        title("Phase Alignment across all electrodes - Hits")
        imagesc(IntanBehaviour.hitTrace(1).time,1:nElectrodes,reshape(PA.HitRewardz,[],size(PA.HitRewardz,3))); colormap(hot);
        ylabel("Electrodes");xlabel("Time (s)");
        xline(0,'-w','Reward','LabelVerticalAlignment','top');caxis([-4 8]);
        xline(-mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Cue Time','LabelVerticalAlignment','top');
        subplot(2,1,2);
        title("Phase Alignment across all electrodes - FA")
        imagesc(IntanBehaviour.missTrace(1).time,1:nElectrodes,reshape(PA.FAz,[],size(PA.FAz,3))); colormap(hot);
        ylabel("Electrodes");xlabel("Time (s)");
        xline(0,'-w','No Reward','LabelVerticalAlignment','top');caxis([-4 8]);
        yline(2.56);
    end
end