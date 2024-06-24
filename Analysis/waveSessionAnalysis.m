%% Plotting rasters 
%  Waves during hits vs miss
plotWaveRaster(Waves.wavesHit,Waves.wavesMiss,IntanBehaviour.cueHitTrace,IntanBehaviour.cueMissTrace,parameters);
sgtitle('Wave Rasters for Hits vs Misses');
subplot(4,1,1);
RTTraceTime = (IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviour.cueHitTrace,2),'.r');

%  Waves during Hits vs FA
plotWaveRaster(Waves.wavesHitReward,Waves.wavesFA,IntanBehaviour.hitTrace,IntanBehaviour.missTrace,parameters);
sgtitle('Wave Rasters for Hits vs FA');

%  Waves during Hits vs FA
plotWaveRaster(Waves.wavesMIHit,Waves.wavesMIFA,IntanBehaviour.MIHitTrace,IntanBehaviour.MIFATrace,parameters);
sgtitle('Wave Rasters for Hits vs FA');


%% Wave Direction
% Waves Hits vs Miss
W = vertcat(selectWaves(Waves.wavesHit,1,3000).waveStart);figure();
prop = arrayfun(@(s) s.waveDir, selectWaves(Waves.wavesHit,1,3000), 'UniformOutput', false);
ax1 = subplot(2,1,1);
rasterPlotPropColor(W,prop,[],1);
RTTraceTime = (IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
% plot(RTTraceTime,1:size(IntanBehaviour.cueHitTrace,2),'.r');xlim([1 size(W,2)]);
title('Wave Hits')

W = vertcat(selectWaves(Waves.wavesMiss,1,3000).waveStart);
prop = arrayfun(@(s) s.waveDir, selectWaves(Waves.wavesMiss,1,3000), 'UniformOutput', false);
ax2 = subplot(2,1,2);
rasterPlotPropColor(W,prop,[],1);
title('Wave Miss')
linkaxes([ax1,ax2],'x');

% Plotting Wave Rasters colorcoded by Direction
srt = 1500;
stp = 1700;
figure();
subplot(1,2,1);
[dirCombHit,~] = plotWaveDirection(selectWaves(Waves.wavesHit(1:end),srt,stp),36,[]);
title('Hits');
subplot(1,2,2);
[dirCombMiss,~] = plotWaveDirection(selectWaves(Waves.wavesMiss(1:end),srt,stp),36,[]);
title('Miss');
sgtitle('Wave Direction')

[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);
% Print the results.
disp('Wave Direction')
disp('p-value:');
disp(p);


srt = 1000;stp = 1300;
figure();
subplot(2,3,1);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(Waves.wavesMIHit,srt,stp),36,[]);
title('Baseline: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombOpto,~] = plotWaveDirection(selectWaves(Waves.wavesMIFA,srt,stp),36,[]);
title('Opto: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1300;stp = 1700;
subplot(2,3,2);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(Waves.wavesMIHit,srt,stp),36,[]);
title('Baseline: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombOpto,~] = plotWaveDirection(selectWaves(Waves.wavesMIFA,srt,stp),36,[]);
title('Opto: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1700;stp = 3000;
subplot(2,3,3);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(Waves.wavesMIHit,srt,stp),36,[]);
title('Baseline: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombOpto,~] = plotWaveDirection(selectWaves(Waves.wavesMIFA,srt,stp),36,[]);
title('Opto: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction')



srt = 1000;stp = 1500;
figure();
subplot(2,3,1);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(Waves.wavesHit,srt,stp),36,[]);
title('Baseline: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombOpto,~] = plotWaveDirection(selectWaves(Waves.wavesMiss,srt,stp),36,[]);
title('Opto: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1500;stp = 1800;
subplot(2,3,2);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(Waves.wavesHit,srt,stp),36,[]);
title('Baseline: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombOpto,~] = plotWaveDirection(selectWaves(Waves.wavesMiss,srt,stp),36,[]);
title('Opto: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1700;stp = 2200;
subplot(2,3,3);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(Waves.wavesHit,srt,stp),36,[]);
title('Baseline: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombOpto,~] = plotWaveDirection(selectWaves(Waves.wavesMiss,srt,stp),36,[]);
title('Opto: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction')

 

