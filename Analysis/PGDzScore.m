PGDHits = vertcat(Waves.wavesHit.PGD);
PGDMiss = vertcat(Waves.wavesMiss.PGD);
% PGDAll = vertcat(PGDHits,PGDMiss);

figure(); hold on;
shadedErrorBar(IntanBehaviour.cueHitTrace(1).time,PGDHits,{@mean,@(x) std(x)/sqrt(size(x,1))});
shadedErrorBar(IntanBehaviour.cueHitTrace(1).time,PGDMiss,{@mean,@(x) std(x)/sqrt(size(x,1))});

figure,hold on;
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGDHits,1),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGDHits,1)-(std(PGDHits,0,1)/sqrt(size(PGDHits,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGDHits,1)+(std(PGDHits,0,1)/sqrt(size(PGDHits,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);

plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGDMiss,1),100,'sgolay',4),'Color', [0.7 0.7 0.7],'LineWidth',1.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGDMiss,1)-(std(PGDMiss,0,1)/sqrt(size(PGDMiss,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGDMiss,1)+(std(PGDMiss,0,1)/sqrt(size(PGDMiss,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','Misses');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');

PGDHits = vertcat(Waves.wavesMIHit.PGD);
PGDMiss = vertcat(Waves.wavesMIFA.PGD);
% PGDAll = vertcat(PGDHits,PGDMiss);

figure(); hold on;
shadedErrorBar(IntanBehaviour.cueHitTrace(1).time,PGDHits,{@mean,@(x) std(x)/sqrt(size(x,1))});
shadedErrorBar(IntanBehaviour.cueHitTrace(1).time,PGDMiss,{@mean,@(x) std(x)/sqrt(size(x,1))});


%%
n = 10000;
PGDHitsShuffled = zeros(n,size(PGDHits,2));
PGDMissShuffled = zeros(n,size(PGDMiss,2));

nHits = size(PGDHits,1);
nMiss = size(PGDMiss,1);

for i=1:n
    r = randperm(nHits+nMiss);
    PGDHitsShuffled(i,:) = mean(PGDAll(r(1:nHits),:),1);
    PGDMissShuffled(i,:) = mean(PGDAll(r(nHits+1:nHits+nMiss),:),1);
end

uHitz = mean(PGDHitsShuffled,1);
stdHitz = std(PGDHitsShuffled,0,1);
PGDHitsz = (mean(PGDHits,1)-uHitz)./stdHitz;

uMissz = mean(PGDMissShuffled,1);
stdMissz = std(PGDMissShuffled,0,1);
PGDMissz = (mean(PGDMiss,1)-uMissz)./stdMissz;


figure(); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,smooth(PGDHitsz,50,'sgolay',4),'-r','LineWidth',1.2); hold on;
% plot(IntanBehaviour.cueMissTrace(1).time,smooth(PGDMissz,50,'sgolay',5),'-k','LineWidth',1);
% plot(IntanBehaviour.cueHitTrace(1).time,PGD.avgPGDFA,'-b','LineWidth',1);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','Misses');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');
