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

figure,hold on;
vq = interp1(interval:interval:interval*nPoints,waveNetDir(1,:),50:50:3000,'spline');
h1 = plot(50:50:3000,vq,'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
vq = interp1(interval:interval:interval*nPoints,waveNetDir(3,:),50:50:3000,'spline');
h2 = plot(50:50:3000,vq,'Color', [0.7 0.7 0.7],'LineWidth',1.5);
vq = interp1(interval:interval:interval*nPoints,waveNetDir(5,:),50:50:3000,'spline');
h3 = plot(50:50:3000,vq,'Color', [0 0 1 0.4],'LineWidth',1.5);
vq = interp1(interval:interval:interval*nPoints,waveNetDir(7,:),50:50:3000,'spline');
h4 = plot(50:50:3000,vq,'Color', [1 0 0 0.4],'LineWidth',1.5);
xline(1501,'--r','Cue');xlabel('Time (ms)'); ylabel('Net Wave Directionality');
legend([h1 h2 h3 h4],'Hits','Miss','MIFAs','MIHits','Location','best'); ylim([0 0.6]);
xlim([1000 3000]);set(gca,'fontsize',14,'linewidth',1.5)

