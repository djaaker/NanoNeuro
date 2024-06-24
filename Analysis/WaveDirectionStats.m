%% Wave direction stats - Hits vs Misses
dirCombHit = horzcat(WavesBaseline.wavesHit(1:end).waveDir);
dirCombMiss = horzcat(WavesBaseline.wavesMiss(1:end).waveDir);

[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);
% Print the results.
disp('Wave Direction')
disp('p-value:');
disp(p);

% Fit VonMisesMixture Model 
nComponents = 2;    
fittedVmBaseline = fitmvmdist(dirCombHit',nComponents,'MaxIter', 250);
fittedVmOpto = fitmvmdist(dirCombMiss',nComponents,'MaxIter', 250);
angles = linspace(-pi, pi, 1000)';
likelihoodsBaseline = fittedVmBaseline.pdf(angles);
likelihoodsOpto = fittedVmOpto.pdf(angles);

figure('Name','Polar Histogram for wave direction');
subplot(1,2,1);
h=polarhistogram(dirCombHit,18,'Normalization','pdf','FaceColor',[153/255 153/255 255/255],'FaceAlpha',1,'EdgeAlpha',0.3);hold on;
polarplot(angles,likelihoodsBaseline,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);rlim([0 0.4]);
polarplot([fittedVmBaseline.mu';fittedVmBaseline.mu'],[0 0;fittedVmBaseline.pdf(fittedVmBaseline.mu)'],'-r','LineWidth',1.5);
title('Wave Direction : Hit');box off;
subplot(1,2,2);
polarhistogram(dirCombMiss,18,'Normalization','pdf','FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.3,'EdgeAlpha',0.3);hold on;
polarplot(angles,likelihoodsOpto,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);rlim([0 0.4]);
polarplot([fittedVmOpto.mu';fittedVmOpto.mu'],[0 0;fittedVmOpto.pdf(fittedVmOpto.mu)'],'-r','LineWidth',1.5);
title('Wave Direction : Miss');box off;

%% Wave direction stats - Hits vs FA
dirCombHit = horzcat(WavesBaseline.wavesHit(1:end).waveDir);
dirCombFA = horzcat(WavesBaseline.wavesFA(1:end).waveDir);

[p,~,~] = circ_kuipertest(dirCombHit, dirCombFA,60,0);
% Print the results.
disp('Wave Direction')
disp('p-value:');
disp(p);

% Fit VonMisesMixture Model 
nComponents = 2;    
fittedVmBaseline = fitmvmdist(dirCombHit',nComponents,'MaxIter', 250);
fittedVmOpto = fitmvmdist(dirCombFA',nComponents,'MaxIter', 250);
angles = linspace(-pi, pi, 1000)';
likelihoodsBaseline = fittedVmBaseline.pdf(angles);
likelihoodsOpto = fittedVmOpto.pdf(angles);

figure('Name','Polar Histogram for wave direction');
subplot(1,2,1);
h=polarhistogram(dirCombHit,18,'Normalization','pdf','FaceColor',[153/255 153/255 255/255],'FaceAlpha',1,'EdgeAlpha',0.3);hold on;
polarplot(angles,likelihoodsBaseline,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);rlim([0 0.4]);
polarplot([fittedVmBaseline.mu';fittedVmBaseline.mu'],[0 0;fittedVmBaseline.pdf(fittedVmBaseline.mu)'],'-r','LineWidth',1.5);
title('Wave Direction : Hit');box off;
subplot(1,2,2);
polarhistogram(dirCombFA,18,'Normalization','pdf','FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.3,'EdgeAlpha',0.3);hold on;
polarplot(angles,likelihoodsOpto,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);rlim([0 0.4]);
polarplot([fittedVmOpto.mu';fittedVmOpto.mu'],[0 0;fittedVmOpto.pdf(fittedVmOpto.mu)'],'-r','LineWidth',1.5);
title('Wave Direction : FA');box off;