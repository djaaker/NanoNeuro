Spikes = SpikesBaseline;
IntanBehaviour = IntanBehaviourBaseline;
LFP = LFPBaseline;
%% Removing clusters that are not correct
Spikes = rejectSpikes(Spikes,0.75,1,parameters);

%% Plotting spikes raster and FR to check 
figure,hold on
for n = 33
subplot(2,1,[1]),Show_Spikes(Spikes.PSTH.hit.spks{n}),axis off
subplot(2,1,[2]),bar(-1500:1500,smoothdata(Spikes.PSTH.hit.spkRates(n,:)),'FaceColor',[28/255 117/255 188/255],'EdgeColor','none')
axis tight, box off, set(gca,'TickDir','out');sgtitle(['Neuron - ' num2str(n)])
set(gca,'fontsize',16)
xline(0)
% xline(mean(IntanBehaviour.reactionTime)*1000)
end
%% Segmenting spiking data according to behaviour
FsScaling = (Spikes.Clusters(1).cluster(1)/Spikes.Clusters(1).spikeTime(1))/parameters.Fs;
for i=1:size(Spikes.Clusters,2)
    Spikes.Clusters(i).spikeIndexLFP = round(Spikes.Clusters(i).cluster/FsScaling);
end
% Making spikes rasters 
Spikes.spikeRaster = zeros(size(Spikes.Clusters,2),size(IntanBehaviour.time,2));
for i=1:size(Spikes.Clusters,2)
    Spikes.spikeRaster(i,Spikes.Clusters(i).spikeIndexLFP) = 1;
end

% Adding spike raster to IntanBehaviour
for i = 1:size(IntanBehaviour.cueHitTrace,2)
    IntanBehaviour.cueHitTrace(i).spikeRaster = Spikes.spikeRaster(:,IntanBehaviour.cueHitTrace(i).LFPIndex(1):IntanBehaviour.cueHitTrace(i).LFPIndex(end));
end
for i = 1:size(IntanBehaviour.cueMissTrace,2)
    IntanBehaviour.cueMissTrace(i).spikeRaster = Spikes.spikeRaster(:,IntanBehaviour.cueMissTrace(i).LFPIndex(1):IntanBehaviour.cueMissTrace(i).LFPIndex(end));
end
for i = 1:size(IntanBehaviour.hitTrace,2)
    IntanBehaviour.hitTrace(i).spikeRaster = Spikes.spikeRaster(:,IntanBehaviour.hitTrace(i).LFPIndex(1):IntanBehaviour.hitTrace(i).LFPIndex(end));
end
for i = 1:size(IntanBehaviour.missTrace,2)
    IntanBehaviour.missTrace(i).spikeRaster = Spikes.spikeRaster(:,IntanBehaviour.missTrace(i).LFPIndex(1):IntanBehaviour.missTrace(i).LFPIndex(end));
end
for i = 1:size(IntanBehaviour.MIHitTrace,2)
    IntanBehaviour.MIHitTrace(i).spikeRaster = Spikes.spikeRaster(:,IntanBehaviour.MIHitTrace(i).LFPIndex(1):IntanBehaviour.MIHitTrace(i).LFPIndex(end));
end
for i = 1:size(IntanBehaviour.MIFATrace,2)
    IntanBehaviour.MIFATrace(i).spikeRaster = Spikes.spikeRaster(:,IntanBehaviour.MIFATrace(i).LFPIndex(1):IntanBehaviour.MIFATrace(i).LFPIndex(end));
end

%% Getting average phase map for each spike
for i=1:Spikes.nSpikes
    Spikes.spikeTrigPhase(i).rawPhaseMap = angle(LFP.xgp(:,:,Spikes.Clusters(i).spikeIndexLFP));
    Spikes.spikeTrigPhase(i).rawAvgPhaseMap = circ_mean(angle(LFP.xgp(:,:,Spikes.Clusters(i).spikeIndexLFP)),[],3);
    Spikes.spikeTrigPhase(i).PA = abs(circ_r(angle(LFP.xgp(:,:,Spikes.Clusters(i).spikeIndexLFP)),[],[],3));
    Spikes.spikeTrigPhase(i).meanPA = mean(Spikes.spikeTrigPhase(i).PA,'all','omitnan');
    a = 1*exp(1i*Spikes.spikeTrigPhase(i).rawAvgPhaseMap);
    Spikes.spikeTrigPhase(i).phaseMap = angle(inpaint_nans(a,3));
    [Spikes.spikeTrigPhase(i).pm,Spikes.spikeTrigPhase(i).pd,Spikes.spikeTrigPhase(i).dx,Spikes.spikeTrigPhase(i).dy] = getPhaseGradient( Spikes.spikeTrigPhase(i).phaseMap, parameters.xspacing,parameters.yspacing);
    Spikes.spikeTrigPhase(i).sourcePoint = find_source_points2(1,parameters.X,parameters.Y,Spikes.spikeTrigPhase(i).dx,Spikes.spikeTrigPhase(i).dy);
    Spikes.spikeTrigPhase(i).PGD = phase_gradient_directionality(Spikes.spikeTrigPhase(i).pm,Spikes.spikeTrigPhase(i).dx,Spikes.spikeTrigPhase(i).dy);
    [Spikes.spikeTrigPhase(i).rho,~,~] = phase_correlation_distance( Spikes.spikeTrigPhase(i).phaseMap,Spikes.spikeTrigPhase(i).sourcePoint, parameters.xspacing,parameters.yspacing );
    [Spikes.spikeTrigPhase(i).vx, Spikes.spikeTrigPhase(i).vy] = wavefront_direction(Spikes.spikeTrigPhase(i).pd,[]);
    Spikes.spikeTrigPhase(i).waveDir = atan2(Spikes.spikeTrigPhase(i).vy,Spikes.spikeTrigPhase(i).vx);
end
%% Plotting spike triggered phase maps
figure,
n = 42;
imagesc(Spikes.spikeTrigPhase(n).phaseMap);
% map = colorcet( 'C2' );
% map = circshift(map,1);
% colormap(map)
colormap("bone");
c = colorbar;
hold on;
[XX,YY] = meshgrid( 1:size(Spikes.spikeTrigPhase(n).dx+1i*Spikes.spikeTrigPhase(n).dy,2), 1:size(Spikes.spikeTrigPhase(n).dx+1i*Spikes.spikeTrigPhase(n).dy,1) );
M = real( exp( 1i * angle(Spikes.spikeTrigPhase(n).dx+1i*Spikes.spikeTrigPhase(n).dy) ) ); N = imag( exp( 1i * angle(Spikes.spikeTrigPhase(n).dx+1i*Spikes.spikeTrigPhase(n).dy) ) );
quiver( XX, YY, M, N, 0.25, 'r' );
scatter(Spikes.spikeTrigPhase(n).sourcePoint(1),Spikes.spikeTrigPhase(n).sourcePoint(2),'filled');

% Plotting PA/SPI 
allPA = cell2mat(arrayfun(@(s) reshape(s.PA,parameters.rows*parameters.cols,[]),Spikes.spikeTrigPhase,'UniformOutput',false));
figure;
imagesc(removeNaNRows(allPA));
colormap hot
c = colorbar;
c.Label.String = 'SPI';
ylabel('Grid Channel')
xlabel('Spiking #')
set(gca,'fontsize',14,'linewidth',1.5)

% Plotting preferred phase 
allPhase = cell2mat(arrayfun(@(s) reshape(s.phaseMap,parameters.rows*parameters.cols,[]),Spikes.spikeTrigPhase,'UniformOutput',false));
figure;
imagesc(removeNaNRows(allPhase));
map = colorcet( 'C2' );map = circshift(map,1);colormap(map)
c = colorbar;
c.Label.String = 'Phase';
ylabel('Grid Channel')
xlabel('Spiking #');title('Preferred Phase');
set(gca,'fontsize',14,'linewidth',1.5)

%% Plotting Phase Gradient directions for the spike triggered phase maps
PMGDirection = cell2mat(arrayfun(@(s) s.waveDir, Spikes.spikeTrigPhase,'UniformOutput',false));
figure,plotDirectionHistogram(PMGDirection,36,[]);
phaseBoundary1 = pi/2;
PMGDirectionRotated = angle(exp(1i*PMGDirection)*exp(1i*-phaseBoundary1));
% figure,plotDirectionHistogram(PMGDirectionRotated,36,[]);
for i=1:Spikes.nSpikes
    if (PMGDirectionRotated(i)>0)
        Spikes.spikeTrigPhase(i).directionCluster = 1;
    else
        Spikes.spikeTrigPhase(i).directionCluster = 2;
    end
end

%% Looking at phase maps just for Hits pre and post cue
pre1 = 1200; pre2 = 1500;
post1 = 1501; post2 = 1801;
for i=1:Spikes.nSpikes % iterrating over all spikes 
    % Pre cue
    a = arrayfun(@(s) angle(s.xgp(:,:,find(s.spikeRaster(i,pre1:pre2)))), IntanBehaviour.cueHitTrace,'UniformOutput',false); 
    Spikes.spikeTrigPhaseMapPre(i).allPMsRaw = cat(3,a{1,1:end});
    Spikes.spikeTrigPhaseMapPre(i).allPMs = inpaintNaNsPM(Spikes.spikeTrigPhaseMapPre(i).allPMsRaw);
    Spikes.spikeTrigPhaseMapPre(i).rawPhaseMap = circ_mean(Spikes.spikeTrigPhaseMapPre(i).allPMsRaw,[],3);
    Spikes.spikeTrigPhaseMapPre(i).PA = abs(circ_r(angle(Spikes.spikeTrigPhaseMapPre(i).allPMsRaw),[],[],3));
    Spikes.spikeTrigPhaseMapPre(i).meanPA = mean(Spikes.spikeTrigPhaseMapPre(i).PA,'all','omitnan');
    Spikes.spikeTrigPhaseMapPre(i).phaseMap = circ_mean(Spikes.spikeTrigPhaseMapPre(i).allPMs,[],3);
    [Spikes.spikeTrigPhaseMapPre(i).pm,Spikes.spikeTrigPhaseMapPre(i).pd,Spikes.spikeTrigPhaseMapPre(i).dx,Spikes.spikeTrigPhaseMapPre(i).dy] = getPhaseGradient( Spikes.spikeTrigPhaseMapPre(i).allPMs, parameters.xspacing,parameters.yspacing);
    Spikes.spikeTrigPhaseMapPre(i).sourcePoint = find_source_points2(1:size(Spikes.spikeTrigPhaseMapPre(i).allPMs,3),parameters.X,parameters.Y,Spikes.spikeTrigPhaseMapPre(i).dx,Spikes.spikeTrigPhaseMapPre(i).dy);
    Spikes.spikeTrigPhaseMapPre(i).PGD = phase_gradient_directionality(Spikes.spikeTrigPhaseMapPre(i).pm,Spikes.spikeTrigPhaseMapPre(i).dx,Spikes.spikeTrigPhaseMapPre(i).dy);
%     [Spikes.spikeTrigPhaseMapPre(i).rho,~,~] = phase_correlation_distance( Spikes.spikeTrigPhaseMapPre(i).phaseMap,Spikes.spikeTrigPhaseMapPre(i).sourcePoint, parameters.xspacing,parameters.yspacing );
    [Spikes.spikeTrigPhaseMapPre(i).vx, Spikes.spikeTrigPhaseMapPre(i).vy] = wavefront_direction(Spikes.spikeTrigPhaseMapPre(i).pd,[]);
    Spikes.spikeTrigPhaseMapPre(i).waveDir = atan2(Spikes.spikeTrigPhaseMapPre(i).vy,Spikes.spikeTrigPhaseMapPre(i).vx);
    
    % Post Cue
    a = arrayfun(@(s) angle(s.xgp(:,:,find(s.spikeRaster(i,post1:post2)))), IntanBehaviour.cueHitTrace,'UniformOutput',false); 
    Spikes.spikeTrigPhaseMapPost(i).allPMsRaw = cat(3,a{1,1:end});
    Spikes.spikeTrigPhaseMapPost(i).allPMs = inpaintNaNsPM(Spikes.spikeTrigPhaseMapPost(i).allPMsRaw);
    Spikes.spikeTrigPhaseMapPost(i).rawPhaseMap = circ_mean(Spikes.spikeTrigPhaseMapPost(i).allPMsRaw,[],3);
    Spikes.spikeTrigPhaseMapPost(i).PA = abs(circ_r(angle(Spikes.spikeTrigPhaseMapPost(i).allPMsRaw),[],[],3));
    Spikes.spikeTrigPhaseMapPost(i).meanPA = mean(Spikes.spikeTrigPhaseMapPost(i).PA,'all','omitnan');
    Spikes.spikeTrigPhaseMapPost(i).phaseMap = circ_mean(Spikes.spikeTrigPhaseMapPost(i).allPMs,[],3);
    [Spikes.spikeTrigPhaseMapPost(i).pm,Spikes.spikeTrigPhaseMapPost(i).pd,Spikes.spikeTrigPhaseMapPost(i).dx,Spikes.spikeTrigPhaseMapPost(i).dy] = getPhaseGradient( Spikes.spikeTrigPhaseMapPost(i).allPMs, parameters.xspacing,parameters.yspacing);
    Spikes.spikeTrigPhaseMapPost(i).sourcePoint = find_source_points2(1:size(Spikes.spikeTrigPhaseMapPost(i).allPMs,3),parameters.X,parameters.Y,Spikes.spikeTrigPhaseMapPost(i).dx,Spikes.spikeTrigPhaseMapPost(i).dy);
    Spikes.spikeTrigPhaseMapPost(i).PGD = phase_gradient_directionality(Spikes.spikeTrigPhaseMapPost(i).pm,Spikes.spikeTrigPhaseMapPost(i).dx,Spikes.spikeTrigPhaseMapPost(i).dy);
%     [Spikes.spikeTrigPhaseMapPost(i).rho,~,~] = phase_correlation_distance( Spikes.spikeTrigPhaseMapPost(i).phaseMap,Spikes.spikeTrigPhaseMapPost(i).sourcePoint, parameters.xspacing,parameters.yspacing );
    [Spikes.spikeTrigPhaseMapPost(i).vx, Spikes.spikeTrigPhaseMapPost(i).vy] = wavefront_direction(Spikes.spikeTrigPhaseMapPost(i).pd,[]);
    Spikes.spikeTrigPhaseMapPost(i).waveDir = atan2(Spikes.spikeTrigPhaseMapPost(i).vy,Spikes.spikeTrigPhaseMapPost(i).vx);
end

%% Plotting

% Plotting PA/SPI 
allPAPre = cell2mat(arrayfun(@(s) reshape(s.PA,parameters.rows*parameters.cols,[]),Spikes.spikeTrigPhaseMapPre,'UniformOutput',false));
figure;
subplot(2,1,1)
imagesc(removeNaNRows(allPAPre));
colormap hot;c = colorbar;c.Label.String = 'SPI';
ylabel('LFP Electrode Channel');xlabel('Spiking #');
set(gca,'fontsize',14,'linewidth',1.5);title('Pre-Cue');
allPAPost = cell2mat(arrayfun(@(s) reshape(s.PA,parameters.rows*parameters.cols,[]),Spikes.spikeTrigPhaseMapPost,'UniformOutput',false));
subplot(2,1,2)
imagesc(removeNaNRows(allPAPost));
colormap hot;c = colorbar;c.Label.String = 'SPI';
ylabel('LFP Electrode Channel');xlabel('Spiking #');
set(gca,'fontsize',14,'linewidth',1.5)
title('Post-Cue');
sgtitle('SPI for Pre-cue and Post-cue');

% Plotting preferred phase 
allPhasePre = cell2mat(arrayfun(@(s) reshape(s.phaseMap,parameters.rows*parameters.cols,[]),Spikes.spikeTrigPhaseMapPre,'UniformOutput',false));
figure;
subplot(2,1,1);
imagesc(removeNaNRows(allPhasePre));
map = colorcet( 'C2' );map = circshift(map,1);colormap(map)
c = colorbar;c.Label.String = 'Phase';
ylabel('Grid Channel');xlabel('Spiking #');title('Pre-Cue');
set(gca,'fontsize',14,'linewidth',1.5);
allPhasePost = cell2mat(arrayfun(@(s) reshape(s.phaseMap,parameters.rows*parameters.cols,[]),Spikes.spikeTrigPhaseMapPost,'UniformOutput',false));
subplot(2,1,2);
imagesc(removeNaNRows(allPhasePost));
map = colorcet( 'C2' );map = circshift(map,1);colormap(map)
c = colorbar;c.Label.String = 'Phase';
ylabel('Grid Channel');xlabel('Spiking #');title('Post-Cue');
set(gca,'fontsize',14,'linewidth',1.5);
sgtitle('Preferred Phase for Pre-cue and Post-cue');

% Plotting PMG direction histograms 
PMGDirectionPre = cell2mat(arrayfun(@(s) circ_mean(s.waveDir,[],2), Spikes.spikeTrigPhaseMapPre,'UniformOutput',false));
PMGDirectionPost = cell2mat(arrayfun(@(s) circ_mean(s.waveDir,[],2), Spikes.spikeTrigPhaseMapPost,'UniformOutput',false));
% PMGDirectionPre = cell2mat(arrayfun(@(s) s.waveDir, Spikes.spikeTrigPhaseMapPre,'UniformOutput',false));
% PMGDirectionPost = cell2mat(arrayfun(@(s) s.waveDir, Spikes.spikeTrigPhaseMapPost,'UniformOutput',false));
figure;
subplot(1,2,1);
title('Pre-Cue');
plotDirectionHistogram(PMGDirectionPre,36,[]);
subplot(1,2,2);
title('Pre-Cue');
plotDirectionHistogram(PMGDirectionPost,36,[]);

% phaseBoundary1 = pi/2;
% PMGDirectionRotated = angle(exp(1i*PMGDirection)*exp(1i*-phaseBoundary1));
% figure,plotDirectionHistogram(PMGDirectionRotated,36,[]);
% for i=1:Spikes.nSpikes
%     if (PMGDirectionRotated(i)>0)
%         Spikes.spikeTrigPhase(i).directionCluster = 1;
%     else
%         Spikes.spikeTrigPhase(i).directionCluster = 2;
%     end
% end

%% Plotting PMG direction for a spike 
n = 44;
figure;
subplot(1,2,1)
plotDirectionHistogram(Spikes.spikeTrigPhaseMapPre(n).waveDir,36,[]);
subplot(1,2,2)
plotDirectionHistogram(Spikes.spikeTrigPhaseMapPost(n).waveDir,36,[]);
%% Wave triggered spiking patterns




