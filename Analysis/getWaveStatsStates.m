function [wavesStat] = getWaveStatsStates(Waves,parameters,plot)

rows = parameters.rows;
cols = parameters.cols;

% Number of waves

WavesPerTrialRest = mean( horzcat(Waves.wavesRest(1:end).nWaves),'all');
totalWavesRest =  sum( horzcat(Waves.wavesRest(1:end).nWaves));
WavesPerTrialRun = mean( horzcat(Waves.wavesRun(1:end).nWaves),'all');
totalWavesRun =  sum( horzcat(Waves.wavesRun(1:end).nWaves));
WavesPerTrialInit = mean( horzcat(Waves.wavesInit(1:end).nWaves),'all');
totalWavesInit =  sum( horzcat(Waves.wavesInit(1:end).nWaves));
WavesPerTrialTerm= mean( horzcat(Waves.wavesTerm(1:end).nWaves),'all');
totalWavesTerm =  sum( horzcat(Waves.wavesTerm(1:end).nWaves));


% Wave Speed stats

speedCombRest = horzcat(Waves.wavesRest(1:end).speed);
avgSpeedRest = mean(speedCombRest);
speedCombRun = horzcat(Waves.wavesRun(1:end).speed);
avgSpeedRun = mean(speedCombRun);
speedCombInit = horzcat(Waves.wavesInit(1:end).speed);
avgSpeedInit = mean(speedCombInit);
speedCombTerm = horzcat(Waves.wavesTerm(1:end).speed);
avgSpeedTerm = mean(speedCombTerm);

% % Perform the t-test.
% [p, t] = ranksum(speedComb, speedCombStop);
% % Print the results.
% disp('Wave Speed')
% disp('h-statistic:');
% disp(t);
% disp('p-value:');
% disp(p);


if plot == 1
    figure('Name','Histogram of wave speeds');
    subplot(4,1,1);
    histfit(speedCombRest,100,'kernel');
    xline(avgSpeedRest,'-r',{'Mean speed = ' num2str(avgSpeedRest) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Rest');xlim([0 inf]);
    subplot(4,1,2);
    histfit(speedCombInit,100,'kernel');
    xline(avgSpeedInit,'-r',{'Mean speed = ' num2str(avgSpeedInit) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Initiation');xlim([0 inf]);
    subplot(4,1,3);
    histfit(speedCombRun,100,'kernel');
    xline(avgSpeedRun,'-r',{'Mean speed = ' num2str(avgSpeedRun) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Run');xlim([0 inf]);
    subplot(4,1,4);
    histfit(speedCombTerm,100,'kernel');
    xline(avgSpeedTerm,'-r',{'Mean speed = ' num2str(avgSpeedTerm) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Motion Term');xlim([0 inf]);

    figure('Name','Wave speeds in Rest, Initiation, Running and Termination');
    group = [ones(size(speedCombRest')); 2.*ones(size(speedCombInit')); 3.*ones(size(speedCombRun')); 4.*ones(size(speedCombTerm'))];
    boxplot([speedCombRest';speedCombInit';speedCombRun';speedCombTerm'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'xtick',1:4,'XTickLabel',{'Rest','Initiation','Run','Termination'});
    ylabel('Wave speed in cm/s');
end


% Wavelength stats
lCombRest = horzcat(Waves.wavesRest(1:end).wavelength);
avglRest = mean(lCombRest);
lCombRun = horzcat(Waves.wavesRun(1:end).wavelength);
avglRun = mean(lCombRun);
lCombInit = horzcat(Waves.wavesInit(1:end).wavelength);
avglInit = mean(lCombInit);
lCombTerm = horzcat(Waves.wavesTerm(1:end).wavelength);
avglTerm = mean(lCombTerm);

if plot == 1
    figure('Name','Histogram of Wavelength');
    subplot(4,1,1);
    histfit(lCombRest,100,'kernel');
    xline(avglRest,'-r',{'Mean l = ' num2str(avglRest) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength : Rest');xlim([0 inf]);
    subplot(4,1,2);
    histfit(lCombInit,100,'kernel');
    xline(avglInit,'-r',{'Mean l = ' num2str(avglInit) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength : Initiation');xlim([0 inf]);
    subplot(4,1,3);
    histfit(lCombRun,100,'kernel');
    xline(avglRun,'-r',{'Mean l = ' num2str(avglRun) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength : Run');xlim([0 inf]);
    subplot(4,1,4);
    histfit(lCombTerm,100,'kernel');
    xline(avglTerm,'-r',{'Mean l = ' num2str(avglTerm) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength : Motion Term');xlim([0 inf]);

    figure('Name','Wavelength in Rest, Initiation, Running and Termination');
    group = [ones(size(lCombRest')); 2.*ones(size(lCombInit')); 3.*ones(size(lCombRun')); 4.*ones(size(lCombTerm'))];
    boxplot([lCombRest';lCombInit';lCombRun';lCombTerm'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'xtick',1:4,'XTickLabel',{'Rest','Initiation','Run','Termination'});
    ylabel('Wavelength in cm');
end



% Wave direction stats

dirCombRest = horzcat(Waves.wavesRest(1:end).waveDir);
avgdirRest = mean(dirCombRest);
dirCombRun = horzcat(Waves.wavesRun(1:end).waveDir);
avgdirRun = mean(dirCombRun);
dirCombInit = horzcat(Waves.wavesInit(1:end).waveDir);
avgdirInit = mean(dirCombInit);
dirCombTerm = horzcat(Waves.wavesTerm(1:end).waveDir);
avgdirTerm = mean(dirCombTerm);


if plot == 1
    figure('Name','Histogram of wave direction');
    subplot(2,2,1);
    polarhistogram(dirCombRest,30);
    title('Wave Direction : Rest');
    subplot(2,2,2);
    polarhistogram(dirCombInit,30);
    title('Wave Direction : Motion Initiation');
    subplot(2,2,3);
    polarhistogram(dirCombRun,30);
    title('Wave Direction : Run');
    subplot(2,2,4);
    polarhistogram(dirCombTerm,30);
    title('Wave Direction : Motion Termination');

    figure('Name','Wave direction in Rest, Initiation, Running and Termination');
    group = [ones(size(dirCombRest')); 2.*ones(size(dirCombInit')); 3.*ones(size(dirCombRun')); 4.*ones(size(dirCombTerm'))];
    boxplot([rad2deg(dirCombRest)';rad2deg(dirCombInit)';rad2deg(dirCombRun)';rad2deg(dirCombTerm)'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'xtick',1:4,'XTickLabel',{'Rest','Initiation','Run','Termination'});
    ylabel('Wave direction in degrees');
end

% Wave source points stats
% sourceComb = horzcat(Waves.wavesStart(1:end).source);
% sourceDen = zeros(rows,cols);
% sourceCombStop = horzcat(Waves.wavesStop(1:end).source);
% sourceDenStop = zeros(rows,cols);
% 
% for j=1:size(sourceComb,2)
%     sourceDen(sourceComb(2,j),sourceComb(1,j)) = sourceDen(sourceComb(2,j),sourceComb(1,j)) + 1;
% end
% maxSourcePoint = max(sourceComb);
% 
% for j=1:size(sourceCombStop,2)
%     sourceDenStop(sourceCombStop(2,j),sourceCombStop(1,j)) = sourceDenStop(sourceCombStop(2,j),sourceCombStop(1,j)) + 1;
% end
% maxSourcePointStop = max(sourceCombStop);
% 
% if plot == 1
%     figure('Name','Spatial map of source points in Motion Initiation and Termination'); 
%     subplot(2,1,1);
%     imagesc(sourceDen);set(gca,'YDir','normal');
%     title('Spatial map of sources points : Motion Inititaion'); colorbar;
%     subplot(2,1,2);
%     imagesc(sourceDenStop);set(gca,'YDir','normal');
%     title('Spatial map of sources points : Motion Termination'); colorbar;
% end 

wavesStat.evaluationPointsRest =  horzcat(Waves.wavesRest(1:end).evaluationPoints);
wavesStat.evaluationPointsInit =  horzcat(Waves.wavesInit(1:end).evaluationPoints);
wavesStat.evaluationPointsRun =  horzcat(Waves.wavesRun(1:end).evaluationPoints);
wavesStat.evaluationPointsTerm =  horzcat(Waves.wavesTerm(1:end).evaluationPoints);

wavesStat.lRest = lCombRest;
wavesStat.avglRest = avglRest;
wavesStat.lInit = lCombInit;
wavesStat.avglInit = avglInit;
wavesStat.lRun = lCombRun;
wavesStat.avglRun = avglRun;
wavesStat.lTerm = lCombTerm;
wavesStat.avglTerm = avglTerm;

wavesStat.speedRest = speedCombRest;
wavesStat.avgSpeedRest = avgSpeedRest;
wavesStat.speedInit = speedCombInit;
wavesStat.avgSpeedInit = avgSpeedInit;
wavesStat.speedRun = speedCombRun;
wavesStat.avgSpeedRun = avgSpeedRun;
wavesStat.speedTerm = speedCombTerm;
wavesStat.avgSpeedTerm = avgSpeedTerm;

wavesStat.dirRest = dirCombRest;
wavesStat.avgdirRest = avgdirRest;
wavesStat.dirInit = dirCombInit;
wavesStat.avgdirInit = avgdirInit;
wavesStat.dirRun = dirCombRun;
wavesStat.avgdirRun = avgdirRun;
wavesStat.dirTerm = dirCombTerm;
wavesStat.avgdirTerm = avgdirTerm;

wavesStat.WavesPerTrialRest = WavesPerTrialRest;
wavesStat.totalWavesRest = totalWavesRest;
wavesStat.WavesPerTrialInit = WavesPerTrialInit;
wavesStat.totalWavesInit = totalWavesInit;
wavesStat.WavesPerTrialRun = WavesPerTrialRun;
wavesStat.totalWavesRun = totalWavesRun;
wavesStat.WavesPerTrialTerm = WavesPerTrialTerm;
wavesStat.totalWavesTerm = totalWavesTerm;

end

