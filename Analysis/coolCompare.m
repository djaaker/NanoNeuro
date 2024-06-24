%% Comaparing Behaviour 
[p,h] = ranksum(IntanBehaviourBaseline.reactionTime,IntanBehaviourCool.reactionTime)
plotBox2(IntanBehaviourBaseline.reactionTime,IntanBehaviourCool.reactionTime);
% xL=xlim;
% yL=ylim;
% text(0.995*xL(2),0.995*yL(2),['p-val = ' num2str(p)],'HorizontalAlignment','right','VerticalAlignment','top')
ylabel('Reaction Time (s)'); title('M2 -> M1 Cool');subtitle(['p-val = ' num2str(p)]);
xtix = {'Baseline','Cool'}; xtixloc = [1 2]; set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);set(gca,'TickDir','out','fontsize',14');
set(gca,'TickDir','out','fontsize',14');

%% Comparing PA 
PABaseline = getPA(IntanBehaviourBaseline,0,1,0,parameters,0);
PACool = getPA(IntanBehaviourCool,0,1,0,parameters,0);

figure();
subplot(2,2,[1,2])
plot(IntanBehaviourBaseline.cueHitTrace(1).time,smooth(squeeze(mean(PABaseline.Hit,[1 2],'omitnan')),50,'sgolay',20),'Color',[166/255 14/255 90/255],'LineWidth',1.5); hold on;
plot(IntanBehaviourCool.cueHitTrace(1).time,smooth(squeeze(mean(PACool.Hit,[1 2],'omitnan')),50,'sgolay',20),'Color',[53/255 189/255 206/255],'LineWidth',1.5);
ylabel("Phase Alignment"); xlabel("Time (s)");
xline(0,'--k','Cue','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xline((PABaseline.PAPeakHit/parameters.Fs-parameters.windowBeforeCue),'--b','Peak','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xlim([-0.5 1.5]);box off;legend('Baseline','Cool');legend('boxoff');set(gca,'TickDir','out','fontsize',14');
title("Phase Alignment : M2 -> M1 Cool");set(gca,'TickDir','out','fontsize',14');

subplot(2,2,3)
PABaselinePeak = PABaseline.PAPeakHit;
xgp = arrayfun(@(s) s.xgp, IntanBehaviourBaseline.cueHitTrace, 'UniformOutput', false);
PABaselineAngles = rmmissing(reshape(cell2mat(cellfun(@(s) reshape(angle(s(:,:,PABaselinePeak)),parameters.rows*parameters.cols,1), xgp,'UniformOutput',false)),[],1));
histogram(PABaselineAngles,18,'FaceAlpha',0.7,'FaceColor',[166/255 14/255 90/255],'Normalization','probability','EdgeColor','none');hold on;
xgp = arrayfun(@(s) s.xgp, IntanBehaviourCool.cueHitTrace, 'UniformOutput', false);
PACoolAngles = rmmissing(reshape(cell2mat(cellfun(@(s) reshape(angle(s(:,:,PABaselinePeak)),parameters.rows*parameters.cols,1), xgp,'UniformOutput',false)),[],1));
histogram(PACoolAngles,18,'FaceAlpha',0.75,'FaceColor',[53/255 189/255 206/255],'Normalization','probability','EdgeColor','none');box off;
xlabel('Peak PA Angle');ylabel('Probablitity');title('PA Angle at Peak');set(gca,'TickDir','out','fontsize',14');
subplot(2,2,4)
polarhistogram(PABaselineAngles,18,'FaceAlpha',0.7,'FaceColor',[166/255 14/255 90/255],'Normalization','probability','EdgeColor','none');hold on;
polarhistogram(PACoolAngles,18,'FaceAlpha',0.75,'FaceColor',[53/255 189/255 206/255],'Normalization','probability','EdgeColor','none');box off;
title('PA Angle at Peak');set(gca,'TickDir','out','fontsize',14');

[p,~,~] = circ_kuipertest(PABaselineAngles, PACoolAngles,60,0);
disp(['Peak PA angle p-val = ' num2str(p)]);

% [p,h] = getContStat(squeeze(reshape(PABaseline.Hit,32,1,[])),squeeze(reshape(PACool.Hit,32,1,[])));


%% Comparing PDG 
PGD.avgPGDBaseline = mean(vertcat(WavesBaseline.wavesHit.PGD),1);
PGD.avgPGDCool = mean(vertcat(WavesCool.wavesHit.PGD),1);

figure(); hold on;
plot(IntanBehaviourBaseline.cueHitTrace(1).time,smooth(PGD.avgPGDBaseline,50,'sgolay',20),'Color',[166/255 14/255 90/255],'LineWidth',1.5);
plot(IntanBehaviourCool.cueHitTrace(1).time,smooth(PGD.avgPGDCool,50,'sgolay',20),'Color',[53/255 189/255 206/255],'LineWidth',1.5);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off; legend('Baseline','Cool');legend('boxoff');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');

%% Comparing Wave Properties - Speed 
% Plotting Wave Rasters colorcoded by Speed
figure();
Spikes = vertcat(selectWaves(WavesBaseline.wavesHit,1,3000).waveStart);
prop = arrayfun(@(s) s.speed, selectWaves(WavesBaseline.wavesHit,1,3000), 'UniformOutput', false);
prop2 = arrayfun(@(s) s.waveDuration, selectWaves(WavesBaseline.wavesHit,1,3000), 'UniformOutput', false);
subplot(2,1,1);
rasterPlotPropColor(Spikes,prop,prop2,0);
RTTraceTime = (IntanBehaviourBaseline.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviourBaseline.cueHitTrace,2),'.r'); 
title('Wave Hits - Baseline');set(gca,'TickDir','out','fontsize',14');
Spikes = vertcat(selectWaves(WavesCool.wavesHit,1,3000).waveStart);
prop = arrayfun(@(s) s.speed, selectWaves(WavesCool.wavesHit,1,3000), 'UniformOutput', false);
prop2 = arrayfun(@(s) s.waveDuration, selectWaves(WavesCool.wavesHit,1,3000), 'UniformOutput', false);
subplot(2,1,2);
rasterPlotPropColor(Spikes,prop,prop2,0);
title('Waves Hits - Cool');set(gca,'TickDir','out','fontsize',14');
RTTraceTime = (IntanBehaviourCool.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviourCool.cueHitTrace,2),'.r');xlim([1 3001]);
sgtitle('Wave rasters for Hits - Baseline vs Cool');

% Plotting wave speed histograms
srt = 1500;
stp = 2000;
speedCombBaseline = horzcat(selectWaves(WavesBaseline.wavesHit,srt,stp).speed);
speedCombCool = horzcat(selectWaves(WavesCool.wavesHit,srt,stp).speed);

figure('Name','Histogram of wave speeds in Baseline and Cool');
subplot(2,1,1);
h1 = histfit(speedCombBaseline,100,'lognormal');
h1(1).FaceAlpha=0.7; h1(1).FaceColor=[166/255 14/255 90/255]; h1(1).EdgeColor='none';
xline(mean(speedCombBaseline),'-r',{'Mean speed = ' num2str(mean(speedCombBaseline)) ' cm/s'});
xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Baseline');box off;xlim([0 20]);
set(gca,'TickDir','out','fontsize',14');
subplot(2,1,2);
h2 = histfit(speedCombCool,100,'lognormal');
h2(1).FaceAlpha=0.7; h2(1).FaceColor=[53/255 189/255 206/255]; h2(1).EdgeColor='none';
xline(mean(speedCombCool),'-r',{'Mean speed = ' num2str(mean(speedCombCool)) ' cm/s'});
xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Cool');box off;xlim([0 20]);
set(gca,'TickDir','out','fontsize',14');

figure('Name','Wave speeds in Baseline and Cool');
plotBox2(speedCombBaseline,speedCombCool);box off;
set(gca,'XTickLabel',{'Baseline','Cool'});set(gca,'TickDir','out','fontsize',14');
ylabel('Wave speed in cm/s');ylim([0 30]);
% Perform the t-test.
[p, t] = ranksum(speedCombBaseline, speedCombCool);
% Print the results.
disp('Wave Speed')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


% Waves speed vs time 
for i = 1:30
    st = (i-1)*100 + 1;
    sp = (i)*100 + 1;
    dirCombBaseline = horzcat(WavesBaseline.wavesHit(1:end).speed);
    evalPointsHit = horzcat(WavesBaseline.wavesHit(1:end).evaluationPoints);
    WaveComb(i).Baseline = dirCombBaseline(evalPointsHit >=st & evalPointsHit <= sp);

    dirCombCool = horzcat(WavesCool.wavesHit(1:end).speed);
    evalPointsMiss = horzcat(WavesCool.wavesHit(1:end).evaluationPoints);
    WaveComb(i).Cool = dirCombCool(evalPointsMiss >=st & evalPointsMiss <= sp);
end

figure(); hold on;
avgVarBaseline = arrayfun(@(s) mean(s.Baseline,'all'), WaveComb);
maxVarBaseline = max(avgVarBaseline);
avgVarCool = arrayfun(@(s) mean(s.Cool,'all'), WaveComb);
maxVarCool = max(avgVarCool);
a = arrayfun(@(s) ranksum(s.Baseline,s.Cool) , WaveComb);
pVal = (a<0.05);
significanceX = find(pVal)*100;
significanceY = 1.1*max(maxVarBaseline,maxVarCool)*ones(1,length(significanceX));
h1 = plot(100:100:3000,(avgVarBaseline),'Color', [166/255 14/255 90/255],'LineWidth',1.5);
h2 = plot(100:100:3000,(avgVarCool),'Color', [53/255 189/255 206/255],'LineWidth',1.5);
xline(1501,'--r','Cue');xlabel('Time (ms)'); ylabel('Wave speed (cm/s)');
plot(significanceX,significanceY,'r*');
legend([h1 h2],'Baseline','Cool','Location','best');
xlim([1000 3000]);set(gca,'TickDir','out','fontsize',14')
title('M2 -> Th Cool');

%% Comparing Wave Properties - Direction

% Plotting Wave Rasters colorcoded by Direction
srt = 1000;stp = 1500;
figure();
subplot(2,3,1);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesBaseline.wavesHit,srt,stp),36,[]);
title('Baseline: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombCool,~] = plotWaveDirection(selectWaves(WavesCool.wavesHit,srt,stp),36,[]);
title('Cool: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombCool,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1500;stp = 2000;
subplot(2,3,2);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesBaseline.wavesHit,srt,stp),36,[]);
title('Baseline: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombCool,~] = plotWaveDirection(selectWaves(WavesCool.wavesHit,srt,stp),36,[]);
title('Cool: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombCool,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1800;stp = 3000;
subplot(2,3,3);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesBaseline.wavesHit,srt,stp),36,[]);
title('Baseline: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombCool,~] = plotWaveDirection(selectWaves(WavesCool.wavesHit,srt,stp),36,[]);
title('Cool: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombCool,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction')


Spikes = vertcat(selectWaves(WavesBaseline.wavesHit,1,3000).waveStart);figure();
prop = arrayfun(@(s) s.waveDir, selectWaves(WavesBaseline.wavesHit,1,3000), 'UniformOutput', false);
prop2 = arrayfun(@(s) s.waveDuration, selectWaves(WavesBaseline.wavesHit,1,3000), 'UniformOutput', false);
ax1 = subplot(2,1,1);
rasterPlotPropColor(Spikes,prop,[],1);
RTTraceTime = (IntanBehaviourBaseline.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
% plot(RTTraceTime,1:size(IntanBehaviourBaseline.cueHitTrace,2),'.r');xlim([1 size(Spikes,2)]);
title('Wave Hits - Baseline');set(gca,'TickDir','out','fontsize',14');

Spikes = vertcat(selectWaves(WavesCool.wavesHit,1,3000).waveStart);
prop = arrayfun(@(s) s.waveDir, selectWaves(WavesCool.wavesHit,1,3000), 'UniformOutput', false);
prop2 = arrayfun(@(s) s.waveDuration, selectWaves(WavesCool.wavesHit,1,3000), 'UniformOutput', false);
ax2 = subplot(2,1,2);
rasterPlotPropColor(Spikes,prop,[],1);
title('Waves Hits - Cool');set(gca,'TickDir','out','fontsize',14');
RTTraceTime = (IntanBehaviourCool.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
% plot(RTTraceTime,1:size(IntanBehaviourCool.cueHitTrace,2),'.r'); xlim([1 size(Spikes,2)]);
sgtitle('Wave rasters for Hits - Baseline vs Cool');
linkaxes([ax1,ax2],'x');
%% Plotting Wave Rasters

%  Waves during hits 
plotWaveRaster(WavesBaseline.wavesHit,WavesCool.wavesHit,IntanBehaviourBaseline.cueHitTrace,IntanBehaviourCool.cueHitTrace,parameters);
sgtitle('Wave Rasters for Hits: Baseline vs Cool');
subplot(4,1,1);
RTTraceTime = (IntanBehaviourBaseline.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviourBaseline.cueHitTrace,2),'.r');
set(gca,'TickDir','out','fontsize',14');
subplot(4,1,3);
RTTraceTime = (IntanBehaviourCool.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviourCool.cueHitTrace,2),'.r');
set(gca,'TickDir','out','fontsize',14');

% Waves during misses
plotWaveRaster(WavesBaseline.wavesMiss,WavesCool.wavesMiss,IntanBehaviourBaseline.cueMissTrace,IntanBehaviourCool.cueMissTrace,parameters);
set(gca,'TickDir','out','fontsize',14');

% Waves during FA
plotWaveRaster(WavesBaseline.wavesFA,WavesCool.wavesFA,IntanBehaviourBaseline.missTrace,IntanBehaviourCool.missTrace,parameters);
set(gca,'TickDir','out','fontsize',14');

% Waves during MIFA
plotWaveRaster(WavesBaseline.wavesMIFA,WavesCool.wavesMIFA,IntanBehaviourBaseline.MIFATrace,IntanBehaviourCool.MIFATrace,parameters);
set(gca,'TickDir','out','fontsize',14');

%% Direction Analyis - Hits vs FA

% Plotting Wave Rasters colorcoded by Direction
srt = 1000;stp = 1400;
figure();
subplot(2,3,1);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesBaseline.wavesMIHit,srt,stp),36,[]);
title('Baseline MIHit: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombCool,~] = plotWaveDirection(selectWaves(WavesBaseline.wavesMIFA,srt,stp),36,[]);
title('Baseline MIFA: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombCool,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1400;stp = 1700;
subplot(2,3,2);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesBaseline.wavesMIHit,srt,stp),36,[]);
title('Baseline MIHit: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombCool,~] = plotWaveDirection(selectWaves(WavesBaseline.wavesMIFA,srt,stp),36,[]);
title('Baseline MIFA: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombCool,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1700;stp = 2000;
subplot(2,3,3);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesBaseline.wavesMIHit,srt,stp),36,[]);
title('Baseline MIHit: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombCool,~] = plotWaveDirection(selectWaves(WavesBaseline.wavesMIFA,srt,stp),36,[]);
title('Baseline MIFA: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombCool,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction')


% Plotting Wave Rasters colorcoded by Direction
srt = 1000;stp = 1400;
figure();
subplot(2,3,1);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesCool.wavesMIHit,srt,stp),36,[]);
title('Cool MIHit: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombCool,~] = plotWaveDirection(selectWaves(WavesCool.wavesMIFA,srt,stp),36,[]);
title('Cool MIFA: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombCool,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1400;stp = 1700;
subplot(2,3,2);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesCool.wavesMIHit,srt,stp),36,[]);
title('Cool MIHit: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombCool,~] = plotWaveDirection(selectWaves(WavesCool.wavesMIFA,srt,stp),36,[]);
title('Cool MIFA: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombCool,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1700;stp = 2000;
subplot(2,3,3);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesCool.wavesMIHit,srt,stp),36,[]);
title('Cool MIHit: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombCool,~] = plotWaveDirection(selectWaves(WavesCool.wavesMIFA,srt,stp),36,[]);
title('Cool MIFA: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombCool,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction')
