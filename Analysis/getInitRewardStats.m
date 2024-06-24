function [wavesStat] = getInitRewardStats(Waves,parameters,plot)

rows = parameters.rows;
cols = parameters.cols;

cueTime = (parameters.Fs*parameters.windowBeforeCue) + 1;
initWindow = 0.8*parameters.Fs;
rewardWindow = 0.6*parameters.Fs;

allWavesEvalPoints = vertcat(Waves.wavesHit(1:end).evaluationPoints);
allWavesSpeed = horzcat(Waves.wavesHit(1:end).speed);
allWavesDir = horzcat(Waves.wavesHit(1:end).waveDir);
allWavesWaveLength = horzcat(Waves.wavesHit(1:end).wavelength);
allWavesWaveDuration = horzcat(Waves.wavesHit(1:end).waveDuration);
allWavesSource = horzcat(Waves.wavesHit(1:end).source);

initWavesIndx = find(allWavesEvalPoints<cueTime & allWavesEvalPoints>cueTime-initWindow);
initWavesSpeed = allWavesSpeed(initWavesIndx);
initWavesDir = allWavesDir(initWavesIndx);
initWavesWavelength = allWavesWaveLength(initWavesIndx);
initWavesWaveDuration = allWavesWaveDuration(initWavesIndx);
initWavesSource = allWavesSource(:,initWavesIndx);


rewardWavesIndx = find(allWavesEvalPoints>cueTime & allWavesEvalPoints<cueTime+rewardWindow);
rewardWavesSpeed = allWavesSpeed(rewardWavesIndx);
rewardWavesDir = allWavesDir(rewardWavesIndx);
rewardWavesWavelength = allWavesWaveLength(rewardWavesIndx);
rewardWavesWaveDuration = allWavesWaveDuration(rewardWavesIndx);
rewardWavesSource = allWavesSource(:,rewardWavesIndx);

% Number of waves 
nWavesInit = size(initWavesIndx);
nWavesReward = size(rewardWavesIndx);

% Wave speed stats
avgSpeedInit = mean(initWavesSpeed);
avgSpeedReward = mean(rewardWavesSpeed);

[p,t] = ranksum(initWavesSpeed,rewardWavesSpeed);
disp('Wave Speed')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);

if plot == 1
    figure('Name','Histogram of wave speeds in motion Initiation and Reward');
    subplot(2,1,1);
    histfit(initWavesSpeed,100,'kernel');
    xline(avgSpeedInit,'-r',{'Mean speed = ' num2str(avgSpeedInit) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Motion Initiation');box off;
    subplot(2,1,2);
    histfit(rewardWavesSpeed,100,'kernel');
    xline(avgSpeedReward,'-r',{'Mean speed = ' num2str(avgSpeedReward) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Reward');box off;

    figure('Name','Wave speeds in motion Initiation and reward');
    group = [ones(size(initWavesSpeed')); 2.*ones(size(rewardWavesSpeed'))];
    boxplot([initWavesSpeed';rewardWavesSpeed'],group,'BoxStyle','filled','PlotStyle','compact');box off;
    set(gca,'XTickLabel',{'Initiation','Reward'});
    ylabel('Wave speed in cm/s');
end


% Wavelength stats
avgWavelengthInit = mean(initWavesWavelength);
avgWavelengthReward = mean(rewardWavesWavelength);

[p,t] = ranksum(initWavesWavelength,rewardWavesWavelength);
disp('Wave Wavelength')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);

if plot == 1
    figure('Name','Histogram of Wavelengths in motion Initiation and Reward');
    subplot(2,1,1);
    histfit(initWavesWavelength,100,'kernel');
    xline(avgWavelengthInit,'-r',{'Mean Wavelength = ' num2str(avgWavelengthInit) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength : Motion Initiation');box off;
    subplot(2,1,2);
    histfit(rewardWavesWavelength,100,'kernel');
    xline(avgWavelengthReward,'-r',{'Mean Wavelength = ' num2str(avgWavelengthReward) ' cm'});box off;
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wave Wavelength : Reward');

    figure('Name','Wavelengths in motion Initiation and Reward');
    group = [ones(size(initWavesWavelength')); 2.*ones(size(rewardWavesWavelength'))];
    boxplot([initWavesWavelength';rewardWavesWavelength'],group,'BoxStyle','filled','PlotStyle','compact');box off;
    set(gca,'XTickLabel',{'Initiation','Reward'});
    ylabel('Wave Wavelength in cm');
end


% Wave Duration stats
avgWaveDurationInit = mean(initWavesWaveDuration);
avgWaveDurationReward = mean(rewardWavesWaveDuration);

[p,t] = ranksum(initWavesWaveDuration,rewardWavesWaveDuration);
disp('Wave Wave Duration')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);

if plot == 1
    figure('Name','Histogram of Wave Duration in motion Initiation and Reward');
    subplot(2,1,1);
    histfit(initWavesWaveDuration,100,'kernel');
    xline(avgWaveDurationInit,'-r',{'Mean WaveDuration = ' num2str(avgWaveDurationInit) ' ms'});
    xlabel('Wave Duration in ms');ylabel('Frequency');title('Wave Duration : Motion Initiation');box off;
    subplot(2,1,2);
    histfit(rewardWavesWaveDuration,100,'kernel');
    xline(avgWaveDurationReward,'-r',{'Mean WaveDuration = ' num2str(avgWaveDurationReward) ' ms'});
    xlabel('WaveDuration in ms');ylabel('Frequency');title('Wave Duration : Reward');box off;

    figure('Name','Wave Duration in motion Initiation and Reward');
    group = [ones(size(initWavesWaveDuration')); 2.*ones(size(rewardWavesWaveDuration'))];
    boxplot([initWavesWaveDuration';rewardWavesWaveDuration'],group,'BoxStyle','filled','PlotStyle','compact');box off;
    set(gca,'XTickLabel',{'Initiation','Reward'});
    ylabel('Wave Duration in ms');
end


% Wave direction stats

avgDirInit = mean(initWavesDir);
avgDirReward = mean(rewardWavesDir);

[p, t] = ranksum(mapAngle360(rad2deg(initWavesDir)), mapAngle360(rad2deg(rewardWavesDir)));
% Print the results.
disp('Wave Direction')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Polar Histogram for wave direction in Motion Initiation and Reward');
    subplot(2,1,1);
    polarhistogram(initWavesDir,30);box off;
    title('Wave Direction : Motion Initiation');
    subplot(2,1,2);
    polarhistogram(rewardWavesDir,30);box off;
    title('Wave Direction : Reward');

    figure('Name','Histogram for wave direction in Motion Initiation and Reward');
    subplot(2,1,1);
    histogram(mapAngle360(rad2deg(initWavesDir)),72);
    title('Wave Direction : Motion Initiation');box off;
    subplot(2,1,2);
    histogram(mapAngle360(rad2deg(rewardWavesDir)),72);box off;
    title('Wave Direction : Reward');

    figure('Name','Wave Direction in motion Initiation and Reward');
    group = [ones(size(initWavesDir'));2.*ones(size(rewardWavesDir'))];
    boxplot([mapAngle360(rad2deg(initWavesDir))';mapAngle360(rad2deg(rewardWavesDir))'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'XTickLabel',{'Initiation','Reward'});box off;
end



wavesStat.initWavesIndx = initWavesIndx;
wavesStat.initWavesSpeed = initWavesSpeed;
wavesStat.initWavesDir = initWavesDir;
wavesStat.initWavesWavelength = initWavesWavelength;
wavesStat.initWavesSource = initWavesSource;

wavesStat.rewardWavesIndx = rewardWavesIndx;
wavesStat.rewardWavesSpeed = rewardWavesSpeed;
wavesStat.rewardWavesDir = rewardWavesDir;
wavesStat.rewardWavesWavelength = rewardWavesWavelength;
wavesStat.rewardWavesSource = rewardWavesSource;

wavesStat.nWavesInit = nWavesInit;
wavesStat.nWavesReward = nWavesReward;

wavesStat.avgSpeedInit = avgSpeedInit;
wavesStat.avgSpeedReward = avgSpeedReward;
wavesStat.avgWavelengthInit = avgWavelengthInit;
wavesStat.avgWavelengthReward = avgWavelengthReward;
wavesStat.avgDirInit = avgDirInit;
wavesStat.avgDirReward = avgDirReward;

end

