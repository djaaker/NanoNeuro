function [wavesStat] = getWaveStats(Waves1,Waves2,parameters,plot)

rows = parameters.rows;
cols = parameters.cols;

% Number of waves
WavesPerTrial = mean( horzcat(Waves1(1:end).nWaves),'all');
totalWaves =  sum( horzcat(Waves1(1:end).nWaves));

WavesPerTrialMiss = mean( horzcat(Waves2(1:end).nWaves),'all');
totalWavesMiss =  sum( horzcat(Waves2(1:end).nWaves));

% Wave Speed stats
speedComb = horzcat(Waves1(1:end).speed);
avgSpeed = mean(speedComb);

speedCombMiss = horzcat(Waves2(1:end).speed);
avgSpeedMiss = mean(speedCombMiss);

% Perform the t-test.
[p, t] = ranksum(speedComb, speedCombMiss);
% Print the results.
disp('Wave Speed')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Histogram of wave speeds in Hits and Misses');
    subplot(2,1,1);
    histfit(speedComb,100,'kernel');
    xline(avgSpeed,'-r',{'Mean speed = ' num2str(avgSpeed) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Hits');box off;
    subplot(2,1,2);
    histfit(speedCombMiss,100,'kernel');
    xline(avgSpeedMiss,'-r',{'Mean speed = ' num2str(avgSpeedMiss) ' cm/s'});
    xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Miss');box off;

    figure('Name','Wave speeds in Hits and Misses');
    group = [ones(size(speedComb')); 2.*ones(size(speedCombMiss'))];
    boxplot([speedComb';speedCombMiss'],group,'BoxStyle','filled','PlotStyle','compact');box off;
    set(gca,'XTickLabel',{'Hits','Misses'});
    ylabel('Wave speed in cm/s');
end


% Wavelength stats
lComb = horzcat(Waves1(1:end).wavelength);
avgl = mean(lComb);

lCombMiss = horzcat(Waves2(1:end).wavelength);
avglMiss = mean(lCombMiss);

% Perform the t-test.
[p, t] = ranksum(lComb, lCombMiss);
% Print the results.
disp('Wavelength')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Histogram of wavelength in Hits and Misses');
    subplot(2,1,1);
    histfit(lComb,100,'kernel');
    xline(avgl,'-r',{'Mean wavelength = ' num2str(avgl) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength:   Hits');box off;
    subplot(2,1,2);
    histfit(lCombMiss,100,'kernel');
    xline(avglMiss,'-r',{'Mean wavelength = ' num2str(avglMiss) ' cm'});
    xlabel('Wavelength in cm');ylabel('Frequency');title('Wavelength :   Misses');box off;

    figure('Name','Wavelength in   Hits and Misses');
    group = [ones(size(lComb')); 2.*ones(size(lCombMiss'))];
    boxplot([lComb';lCombMiss'],group,'BoxStyle','filled','PlotStyle','compact');box off;
    set(gca,'XTickLabel',{'Hits','Misses'});
    xlabel('Wavelength in cm');
end

% Wave Duration stats
tComb = horzcat(Waves1(1:end).waveDuration);
avgt = mean(tComb);

tCombMiss = horzcat(Waves2(1:end).waveDuration);
avgtMiss = mean(tCombMiss);

% Perform the t-test.
[p, t] = ranksum(tComb, tCombMiss);
% Print the results.
disp('Wave Duration')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Histogram of wave duration in Hits and Misses');
    subplot(2,1,1);
    histfit(tComb,100,'kernel');
    xline(avgt,'-r',{'Mean wave duration = ' num2str(avgt) ' ms'});
    xlabel('Wave Duration in ms');ylabel('Frequency');title('Wave Duration:   Hits');box off;
    subplot(2,1,2);
    histfit(tCombMiss,100,'kernel');
    xline(avgtMiss,'-r',{'Mean wave duration = ' num2str(avgtMiss) ' ms'});
    xlabel('Wave Duration in ms');ylabel('Frequency');title('Wave Duration :   Misses');box off;

    figure('Name','Wave duration in   Hits and Misses');
    group = [ones(size(tComb')); 2.*ones(size(tCombMiss'))];
    boxplot([tComb';tCombMiss'],group,'BoxStyle','filled','PlotStyle','compact');box off;
    set(gca,'XTickLabel',{'Hits','Misses'});
    xlabel('Wave Duration in ms');
end

% Wave Amplitude stats
ampComb = horzcat(Waves1(1:end).waveAmp);
avgAmp = mean(ampComb);

ampCombMiss = horzcat(Waves2(1:end).waveAmp);
avgAmpMiss = mean(ampCombMiss);

% Perform the t-test.
[p, t] = ranksum(ampComb, ampCombMiss);
% Print the results.
disp('Wave Amplitude')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Histogram of wave amplitudes in Hits and Misses');
    subplot(2,1,1);
    histfit(ampComb,100,'kernel');
    xline(avgAmp,'-r',{'Mean wave amplitude = ' num2str(avgAmp) ' ms'});
    xlabel('Wave Amplitude in ms');ylabel('Frequency');title('Wave Amplitude: Hits');box off; xlim([0 250]);
    subplot(2,1,2);
    histfit(ampCombMiss,100,'kernel');
    xline(avgAmpMiss,'-r',{'Mean wave amplitude = ' num2str(avgAmpMiss) ' ms'});
    xlabel('Wave Amplitude in ms');ylabel('Frequency');title('Wave Amplitude : Misses');box off;xlim([0 250]);

    figure('Name','Wave amplitude in   Hits and Misses');
    group = [ones(size(ampComb')); 2.*ones(size(ampCombMiss'))];
    boxplot([ampComb';ampCombMiss'],group,'BoxStyle','filled','PlotStyle','compact');box off;
    set(gca,'XTickLabel',{'Hits','Misses'});
    xlabel('Wave Amplitude in \muV');
end

% Wave direction stats
dirComb = horzcat(Waves1(1:end).waveDir);
avgDir = mean(dirComb);
dirCombMiss = horzcat(Waves2(1:end).waveDir);
avgDirMiss = mean(dirCombMiss);

[p,~,~] = circ_kuipertest(dirComb, dirCombMiss,60,0);
% Print the results.
disp('Wave Direction')
disp('p-value:');
disp(p);


if plot == 1
    figure('Name','Polar Histogram for wave direction in Hits and Misses');
    subplot(2,1,1);
    polarhistogram(dirComb,60);
    title('Wave Direction :   Hits');box off;
    subplot(2,1,2);
    polarhistogram(dirCombMiss,60);
    title('Wave Direction :   Misses');box off;

    figure('Name','Wave Direction in   Hits and Misses');
    group = [ones(size(dirComb'));2.*ones(size(dirCombMiss'))];
    boxplot([mapAngle360(rad2deg(dirComb))';mapAngle360(rad2deg(dirCombMiss))'],group,'BoxStyle','filled','PlotStyle','compact');
    set(gca,'XTickLabel',{'Hits','Misses'});box off;
end

figure('Name','Histogram for wave direction in   Hits and Misses');
subplot(2,1,1);
histogram(mapAngle360(rad2deg(dirComb)),60);
title('Wave Direction :  Hits');box off;
subplot(2,1,2);
histogram(mapAngle360(rad2deg(dirCombMiss)),60);
title('Wave Direction :  Misses');box off;


% Wave source points stats
sourceComb = vertcat(Waves1(1:end).source);
sourceDen = zeros(rows,cols);
sourceCombMiss = vertcat(Waves2(1:end).source);
sourceDenMiss = zeros(rows,cols);

for j=1:size(sourceComb,2)
    sourceDen(sourceComb(2,j),sourceComb(1,j)) = sourceDen(sourceComb(2,j),sourceComb(1,j)) + 1;
end
maxSourcePoint = max(sourceComb);

for j=1:size(sourceCombMiss,2)
    sourceDenMiss(sourceCombMiss(2,j),sourceCombMiss(1,j)) = sourceDenMiss(sourceCombMiss(2,j),sourceCombMiss(1,j)) + 1;
end
maxSourcePointMiss = max(sourceCombMiss);

if plot == 1
    figure('Name','Spatial map of source points in   Hits and Misses'); 
    subplot(2,1,1);
    imagesc(sourceDen);set(gca,'YDir','normal');box off;
    title('Spatial map of sources points :   Inititaion'); colorbar;
    subplot(2,1,2);
    imagesc(sourceDenMiss);set(gca,'YDir','normal');box off;
    title('Spatial map of sources points :   Misses'); colorbar;
end 

wavesStat.evaluationPoints =  horzcat(Waves1(1:end).evaluationPoints);
wavesStat.evaluationPointsMiss =  horzcat(Waves2(1:end).evaluationPoints);
wavesStat.speed = speedComb;
wavesStat.avgSpeed = avgSpeed;
wavesStat.speedMiss = speedCombMiss;
wavesStat.avgSpeedMiss = avgSpeedMiss;

wavesStat.velDir = dirComb;
wavesStat.avgDir = avgDir;
wavesStat.velDirMiss = dirCombMiss;
wavesStat.avgDirMiss = avgDirMiss;

wavesStat.sourcePoints = sourceComb;
wavesStat.maxSourcePoint = maxSourcePoint;
wavesStat.sourcePointsMiss = sourceCombMiss;
wavesStat.maxSourcePointMiss = maxSourcePointMiss;

wavesStat.WavesPerTrial = WavesPerTrial;
wavesStat.totalWaves = totalWaves;
wavesStat.WavesPerTrialMiss = WavesPerTrialMiss;
wavesStat.totalWavesMiss = totalWavesMiss;


% Getting wave distribution 



end

