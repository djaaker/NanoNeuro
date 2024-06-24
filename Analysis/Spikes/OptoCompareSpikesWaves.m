
%% Plotting PA/SPI 
allPABaseline = cell2mat(arrayfun(@(s) reshape(s.PA,parameters.rows*parameters.cols,[]),SpikesBaseline.spikeTrigPhase,'UniformOutput',false));
allPAOpto = cell2mat(arrayfun(@(s) reshape(s.PA,parameters.rows*parameters.cols,[]),SpikesOpto.spikeTrigPhase,'UniformOutput',false));
figure;
subplot(2,2,1)
imagesc(removeNaNRows(allPABaseline));
colormap hot;c = colorbar;c.Label.String = 'SPI';
ylabel('Grid Channel');xlabel('Neurons');title('Baseline')
set(gca,'fontsize',14,'linewidth',1.5);box off;
subplot(2,2,2)
imagesc(removeNaNRows(allPAOpto));
colormap hot;c = colorbar;c.Label.String = 'SPI';
ylabel('Grid Channel');xlabel('Neurons');title('Opto')
set(gca,'fontsize',14,'linewidth',1.5);box off;

subplot(2,2,[3,4])
depthBaseline = cell2mat(arrayfun(@(s) s.spikeDepth, SpikesBaseline.Clusters,'UniformOutput',false));
[~,sortspikeDepth1] = sort(depthBaseline);
depthOpto = cell2mat(arrayfun(@(s) s.spikeDepth, SpikesOpto.Clusters,'UniformOutput',false));
[~,sortspikeDepth2] = sort(depthOpto);
plot(depthBaseline(sortspikeDepth1),mean(allPABaseline(:,sortspikeDepth1),1,'omitnan'),'Color',[0.7 0.7 0.7],'LineWidth',1.5); hold on;
plot(depthOpto(sortspikeDepth2),mean(allPAOpto(:,sortspikeDepth2),1,'omitnan'),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5); hold on;
ylabel("SPI"); xlabel('Depth (um)');
xline(0,'--k','Cue','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
box off;legend('Baseline','Opto');legend('boxoff');set(gca,'TickDir','out','fontsize',14');
title("SPI : M2 -> M1 Opto");set(gca,'TickDir','out','fontsize',14');



%% Plotting preferred phase 
allPhaseBaseline = cell2mat(arrayfun(@(s) reshape(s.phaseMap,parameters.rows*parameters.cols,[]),SpikesBaseline.spikeTrigPhase,'UniformOutput',false));
allPhaseOpto = cell2mat(arrayfun(@(s) reshape(s.phaseMap,parameters.rows*parameters.cols,[]),SpikesOpto.spikeTrigPhase,'UniformOutput',false));
figure;
subplot(2,1,1);
imagesc(removeNaNRows(allPhaseBaseline));
map = colorcet( 'C2' );map = circshift(map,1);colormap(map)
c = colorbar;
c.Label.String = 'Phase';
ylabel('Grid Channel')
xlabel('Spiking #');title('Preferred Phase');
set(gca,'fontsize',14,'linewidth',1.5)
subplot(2,1,2);
imagesc(removeNaNRows(allPhaseOpto));
map = colorcet( 'C2' );map = circshift(map,1);colormap(map)
c = colorbar;
c.Label.String = 'Phase';
ylabel('Grid Channel')
xlabel('Spiking #');title('Preferred Phase');
set(gca,'fontsize',14,'linewidth',1.5)

%% Plotting Phase Map Direction
% PMGDirectionBaseline = cell2mat(arrayfun(@(s) circ_mean(s.waveDir,[],2), SpikesBaseline.spikeTrigPhase,'UniformOutput',false));
% PMGDirectionOpto = cell2mat(arrayfun(@(s) circ_mean(s.waveDir,[],2), SpikesOpto.spikeTrigPhase,'UniformOutput',false));
PMGDirectionBaseline = cell2mat(arrayfun(@(s) s.waveDir, SpikesBaseline.spikeTrigPhase,'UniformOutput',false));
PMGDirectionOpto = cell2mat(arrayfun(@(s) s.waveDir, SpikesOpto.spikeTrigPhase,'UniformOutput',false));
figure;
subplot(1,2,1);
title('Baseline');
plotDirectionHistogram(PMGDirectionBaseline,36,[]);
subplot(1,2,2);
title('Opto');
plotDirectionHistogram(PMGDirectionOpto,36,[]);