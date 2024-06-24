%% Plotting Wave figures for a single session 
% Hits vs Miss 
plotWaveRaster(Waves.wavesHit,Waves.wavesMiss,IntanBehaviour.cueHitTrace,IntanBehaviour.cueMissTrace,parameters);
sgtitle('Wave Rasters for Hits vs Misses');
subplot(4,1,1);
RTTraceTime = (IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviour.cueHitTrace,2),'.r');
set(gca,'TickDir','out','fontsize',14');

% Hits vs FA - Reward Aligned 
plotWaveRaster(Waves.wavesHitReward,Waves.wavesFA,IntanBehaviour.hitTrace,IntanBehaviour.missTrace,parameters);
sgtitle('Wave Rasters for Hits vs FAs - Reward aligned');
subplot(4,1,1);
RTTraceTime = -(IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforePull*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviour.hitTrace,2),'.r');
set(gca,'TickDir','out','fontsize',14');

% Hits vs FA - MI Aligned 
plotWaveRaster(Waves.wavesMIHit,Waves.wavesMIFA,IntanBehaviour.MIHitTrace,IntanBehaviour.MIFATrace,parameters);
sgtitle('Wave Rasters for Hits vs FAs - MI aligned');
subplot(4,1,1);
RTTraceTime = -(IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforeMI*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviour.MIHitTrace,2),'.r');
set(gca,'TickDir','out','fontsize',14');

%% Plotting wave direction progression 
% Plotting Wave Rasters colorcoded by Direction - Hits vs Miss
srt = 1;stp = 1500;
figure();
subplot(2,3,1);
[dirCombHit,~] = plotWaveDirection(selectWaves(Waves.wavesHit,srt,stp),36,[]);
title('Hits: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombMiss,~] = plotWaveDirection(selectWaves(Waves.wavesMiss,srt,stp),36,[]);
title('Miss: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1500;stp = 1800;
subplot(2,3,2);
[dirCombHit,~] = plotWaveDirection(selectWaves(Waves.wavesHit,srt,stp),36,[]);
title('Hits: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombMiss,~] = plotWaveDirection(selectWaves(Waves.wavesMiss,srt,stp),36,[]);
title('Miss: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1800;stp = 3000;
subplot(2,3,3);
[dirCombHit,~] = plotWaveDirection(selectWaves(Waves.wavesHit,srt,stp),36,[]);
title('Hits: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombMiss,~] = plotWaveDirection(selectWaves(Waves.wavesMiss,srt,stp),36,[]);
title('Miss: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction - Hits vs Miss')


% Plotting Wave Rasters colorcoded by Direction - Hits vs FA - Reward
% aligned
srt = 1;stp = 1500;
figure();
subplot(2,3,1);
[dirCombHit,~] = plotWaveDirection(selectWaves(Waves.wavesHitReward,srt,stp),36,[]);
title('Hits: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombMiss,~] = plotWaveDirection(selectWaves(Waves.wavesFA,srt,stp),36,[]);
title('FAs: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1500;stp = 1700;
subplot(2,3,2);
[dirCombHit,~] = plotWaveDirection(selectWaves(Waves.wavesHitReward,srt,stp),36,[]);
title('Hits: Reward Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombMiss,~] = plotWaveDirection(selectWaves(Waves.wavesFA,srt,stp),36,[]);
title('FAs: Reward Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1700;stp = 3000;
subplot(2,3,3);
[dirCombHit,~] = plotWaveDirection(selectWaves(Waves.wavesHitReward,srt,stp),36,[]);
title('Hits: Post Reward Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombMiss,~] = plotWaveDirection(selectWaves(Waves.wavesFA,srt,stp),36,[]);
title('FAs: Post Reward Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction: Hits vs FA - Reward aligned')


% Plotting Wave Rasters colorcoded by Direction - Hits vs FA - MI
% aligned
srt = 1;stp = 1500;
figure();
subplot(2,3,1);
[dirCombHit,~] = plotWaveDirection(selectWaves(Waves.wavesMIHit,srt,stp),36,[]);
title('Hits: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombMiss,~] = plotWaveDirection(selectWaves(Waves.wavesMIFA,srt,stp),36,[]);
title('FAs: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1500;stp = 1700;
subplot(2,3,2);
[dirCombHit,~] = plotWaveDirection(selectWaves(Waves.wavesMIHit,srt,stp),36,[]);
title('Hits: MI Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombMiss,~] = plotWaveDirection(selectWaves(Waves.wavesMIFA,srt,stp),36,[]);
title('FAs: MI Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1700;stp = 3000;
subplot(2,3,3);
[dirCombHit,~] = plotWaveDirection(selectWaves(Waves.wavesMIHit,srt,stp),36,[]);
title('Hits: Post MI Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombMiss,~] = plotWaveDirection(selectWaves(Waves.wavesMIFA,srt,stp),36,[]);
title('FAs: Post MI Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction: Hits vs FA - MI aligned')

% Raster color coded by wave direction 

wavePropRaster = vertcat(selectWaves(Waves.wavesHit,1,3000).waveStart);figure();
prop = arrayfun(@(s) s.waveDir, selectWaves(Waves.wavesHit,1,3000), 'UniformOutput', false);
prop2 = arrayfun(@(s) s.waveDuration, selectWaves(Waves.wavesHit,1,3000), 'UniformOutput', false);
ax1 = subplot(2,1,1);
rasterPlotPropColor(wavePropRaster,prop,prop2,1);
RTTraceTime = (IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
% plot(RTTraceTime,1:size(IntanBehaviourBaseline.cueHitTrace,2),'.r');xlim([1 size(wavePropRaster,2)]);
title('Wave Hits - Baseline');set(gca,'TickDir','out','fontsize',14');

wavePropRaster = vertcat(selectWaves(Waves.wavesMiss,1,3000).waveStart);
prop = arrayfun(@(s) s.waveDir, selectWaves(Waves.wavesMiss,1,3000), 'UniformOutput', false);
prop2 = arrayfun(@(s) s.waveDuration, selectWaves(Waves.wavesMiss,1,3000), 'UniformOutput', false);
ax2 = subplot(2,1,2);
rasterPlotPropColor(wavePropRaster,prop,prop2,1);
title('Waves Hits - Opto');set(gca,'TickDir','out','fontsize',14');
sgtitle('Wave rasters for Hits vs Misses');
linkaxes([ax1,ax2],'x');

% Plotting wave direction as function of time
nPoints = 15; interval = (parameters.Fs*(parameters.windowAfterCue+parameters.windowBeforeCue))/nPoints;
waveNetDir = zeros(8,nPoints);
for i=1:nPoints
    st = (i-1)*interval + 1;
    sp = (i)*interval + 1;
    waveNetDir(1,i) = circ_r((horzcat(selectWaves(Waves.wavesHit,st,sp).waveDir))');
    waveNetDir(2,i) = circ_mean((horzcat(selectWaves(Waves.wavesHit,st,sp).waveDir))');
    waveNetDir(3,i) = circ_r((horzcat(selectWaves(Waves.wavesMiss,st,sp).waveDir))');
    waveNetDir(4,i) = circ_mean((horzcat(selectWaves(Waves.wavesMiss,st,sp).waveDir))');
    waveNetDir(5,i) = circ_r((horzcat(selectWaves(Waves.wavesMIFA,st,sp).waveDir))');
    waveNetDir(6,i) = circ_mean((horzcat(selectWaves(Waves.wavesMIFA,st,sp).waveDir))');
    waveNetDir(7,i) = circ_r((horzcat(selectWaves(Waves.wavesMIHit,st,sp).waveDir))');
    waveNetDir(8,i) = circ_mean((horzcat(selectWaves(Waves.wavesMIHit,st,sp).waveDir))');
end

figure;hold on;
h1 = plot(interval:interval:interval*nPoints,waveNetDir(1,:),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
h2 = plot(interval:interval:interval*nPoints,waveNetDir(3,:),'Color', [0.7 0.7 0.7],'LineWidth',1.5);
h3 = plot(interval:interval:interval*nPoints,waveNetDir(5,:),'Color', [0 0 1 0.4],'LineWidth',1.5);
h4 = plot(interval:interval:interval*nPoints,waveNetDir(7,:),'Color', [1 0 0 0.4],'LineWidth',1.5);
% h = scatter(interval:interval:interval*nPoints,waveNetDir(1,:),80,waveNetDir(2,:),'filled');
xline(1501,'--r','Cue');xlabel('Time (ms)'); ylabel('Net Wave Directionality');
legend([h1 h2 h3 h4],'Hits','Miss','MIFAs','MIHits','Location','best'); ylim([0 0.6]);
map = colorcet( 'C2' );map = circshift(map,1);colormap(map)
c = colorbar;c.Label.String = 'Phase';

%% Plotting wave rasters color coded with average frequency of the waves

wavePropRaster = vertcat(selectWaves(Waves.wavesHit,1,3000).waveStart);
figure();
prop = arrayfun(@(s) s.waveFreq, selectWaves(Waves.wavesHit,1,3000), 'UniformOutput', false);
% prop2 = arrayfun(@(s) s.waveDuration, selectWaves(Waves.wavesHit,1,3000), 'UniformOutput', false);
% ax1 = subplot(2,1,1);
rasterPlotPropColor(wavePropRaster,prop,[],0);
RTTraceTime = (IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviourBaseline.cueHitTrace,2),'.r');xlim([1 size(wavePropRaster,2)]);
title('Wave Hits - Baseline');set(gca,'TickDir','out','fontsize',14');
c = colorbar;c.Label.String='Frequency (Hz)';ylabel('Trials');xlabel('Time (ms)');

% Plotting wave freq as function of time
nPoints = 15; interval = (parameters.Fs*(parameters.windowAfterCue+parameters.windowBeforeCue))/nPoints;
waveAvgFreq = zeros(4,nPoints);
for i=1:nPoints
    st = (i-1)*interval + 1;
    sp = (i)*interval + 1;
    waveAvgFreq(1,i) = mean(horzcat(selectWaves(Waves.wavesHit,st,sp).waveFreq),'all','omitnan');
    waveAvgFreq(2,i) = mean(horzcat(selectWaves(Waves.wavesMiss,st,sp).waveFreq),'all','omitnan');
    waveAvgFreq(3,i) = mean(horzcat(selectWaves(Waves.wavesMIFA,st,sp).waveFreq),'all','omitnan');
    waveAvgFreq(4,i) = mean(horzcat(selectWaves(Waves.wavesMIHit,st,sp).waveFreq),'all','omitnan');
end

figure;hold on;
h1 = plot(interval:interval:interval*nPoints,waveAvgFreq(1,:),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
h2 = plot(interval:interval:interval*nPoints,waveAvgFreq(2,:),'Color', [0.7 0.7 0.7],'LineWidth',1.5);
h3 = plot(interval:interval:interval*nPoints,waveAvgFreq(3,:),'Color', [0 0 1 0.4],'LineWidth',1.5);
h4 = plot(interval:interval:interval*nPoints,waveAvgFreq(4,:),'Color', [1 0 0 0.4],'LineWidth',1.5);
% h = scatter(interval:interval:interval*nPoints,waveNetDir(1,:),80,waveNetDir(2,:),'filled');
xline(1501,'--r','Cue');xlabel('Time (ms)'); ylabel('Average Wave Frequency (Hz)');
legend([h1 h2 h3 h4],'Hits','Miss','MIFAs','MIHits','Location','best'); ylim([5 15]);

%% Wave Speed 

% Plotting wave freq as function of time
nPoints = 30; interval = (parameters.Fs*(parameters.windowAfterCue+parameters.windowBeforeCue))/nPoints;
waveAvgFreq = zeros(4,nPoints);
for i=1:nPoints
    st = (i-1)*interval + 1;
    sp = (i)*interval + 1;
    WaveSpeed(i).speedHit = horzcat(selectWaves(Waves.wavesHit,st,sp).speed);
    WaveSpeed(i).speedMiss = horzcat(selectWaves(Waves.wavesMiss,st,sp).speed);
    WaveSpeed(i).speedMIHit = horzcat(selectWaves(Waves.wavesMIHit,st,sp).speed);
    WaveSpeed(i).speedMIFA = horzcat(selectWaves(Waves.wavesMIFA,st,sp).speed);
end

t = interval:interval:interval*nPoints;
y = cell2mat(arrayfun(@(s) mean(s.speedHit,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedHit,0,'all','omitnan'),WaveSpeed,'UniformOutput',false));
figure;hold on;
h1 = errorbar(t,y,err,'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
y = cell2mat(arrayfun(@(s) mean(s.speedMiss,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedMiss,0,'all','omitnan'),WaveSpeed,'UniformOutput',false));
h2 = errorbar(t,y,err,'Color', [0.7 0.7 0.7],'LineWidth',1.5);
y = cell2mat(arrayfun(@(s) mean(s.speedMIHit,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedMIHit,0,'all','omitnan'),WaveSpeed,'UniformOutput',false));
h3 = errorbar(t,y,err,'Color', [0 0 1 0.4],'LineWidth',1.5);
y = cell2mat(arrayfun(@(s) mean(s.speedMIFA,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedMIFA,0,'all','omitnan'),WaveSpeed,'UniformOutput',false));
h4 = errorbar(t,y,err,'Color', [1 0 0 0.4],'LineWidth',1.5);

xline(1501,'--r','Cue');xlabel('Time (ms)'); ylabel('Average Wave Frequency (Hz)');
legend([h1 h2 h3 h4],'Hits','Miss','MIFAs','MIHits','Location','best'); %ylim([5 15]);
