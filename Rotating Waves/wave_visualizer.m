%function wave_visualizer(behaviour,trial)

% parameters
parameters.Fs = 1000; % data sampling rate [Hz]
parameters.filter_order = 4; parameters.f = [40 80]; % filter parameters
parameters.start_time = -300; % time point to start analysis
parameters.stop_time = 0; % stop analysis (relative to target onset)
parameters.pixel_spacing = 0.4; % spacing between electrodes [mm]
parameters.lp = 0; % cutoff for negative frequency detection [Hz]
parameters.evaluation_angle = pi; parameters.tol = 0.2; % evaluation points
parameters.windowBeforeCue = 1.5; % in seconds 
parameters.windowAfterCue = 1.5; % in seconds 
parameters.plot_rho_value = 0.5; plot_time = 20; pause_length = 0.2; % plot settings
parameters.verbose = true; parameters.traces = {behaviour.cueHitTrace, behaviour.cueMissTrace};
parameters.extract = true;
parameters.plot = true;

% datasets = {'Day14_BaselineWaves2.mat'};


if parameters.extract == true 
    rotWaves = struct();
%     for ii = 1:numel(datasets)
%         data = datasets{ii};
%         load(data,'IntanBehaviourBaseline');
%         behaviour = IntanBehaviourBaseline;
        for jj = 1:numel(parameters.traces)
            cueTrace = parameters.traces{jj};
            if jj == 1
                fieldName = 'wavesHit';
            elseif jj == 2
                fieldName = 'wavesMiss';
            end
            condWaves = processWaveTrials(cueTrace, parameters);
            rotWaves.(fieldName) = condWaves;
        end
%     end
end


% rotWaves.wavesHit = wavesHit;

% end  

save('RotWavesHits_4080_D14Cool.mat',"rotWaves")
% load('RotWavesHits.mat','rotWaves')

% rotWaves = load('RotWavesHits_5-40');
hitsVel = []; hitsStart = []; hitsQuad = []; wavesPresAll = NaN(3001,143);
for i = 1:3
    % wavesPresAll(:,i) = rotWaves.wavesHit(i).wavePresent;
    hitsVel = [hitsVel rotWaves.wavesHit(i).angVel];
    hitsStart = [hitsStart rotWaves.wavesHit(i).start];
    hitsQuad = [hitsQuad rotWaves.wavesHit(i).quadrant];
end

if parameters.plot == true
    plotWaveRaster(rotWaves.wavesHit,[],behaviour.cueHitTrace,[])
    
    figure;
    
    % First subplot: scatterplot of angular velocity for hit trials
    subplot(2,1,1);  % Create subplot (2 rows, 1 column, 1st plot)
    scatter(hitsStart, abs(hitsVel), 'filled');
    title('Angular Velocity for Hit Trials');
    xlabel('Time (ms)');
    grid on;
    
    % Remove y-axis for the scatter plot
    set(gca, 'ytick', [])
    
    % Add trend line to scatter plot
    p = polyfit(hitsStart, abs(hitsVel), 1);
    yfit = polyval(p, hitsStart);
    hold on;
    plot(hitsStart, yfit, '-b', 'LineWidth', 1.5);


%     % Add average
%     avgVel = zeros(1,30);
%     for i=1:30
%         pos = find(hitsStart > (i-1)*100 + 1 & hitsStart <= i*100);
%         avgVel(i) =  mean(hitsVel(pos),'all');
%     end
%     hold on;
%     plot(avgVel)
%     hold off;
    
    % Add a red vertical line in the center
    xCenter = (max(hitsStart) + min(hitsStart)) / 2;
    line([xCenter xCenter], get(gca, 'ylim'), 'Color', 'r', 'LineWidth', 1.5);
    
    % Distribution of waves in hit trials
    subplot(2, 1, 2);
    histogram(hitsQuad, 'BinLimits', [1 17], 'BinWidth', 1);
    hold on;
    xlim([0 17]);
    
    % Fit a normal distribution to the data
    pd = fitdist(hitsQuad', 'Normal');
    xValues = linspace(0, 17, 100); % X values from 0 to 16
    yValues = pdf(pd, xValues) * numel(hitsQuad) * (xValues(2) - xValues(1)); % Scaled PDF
    
    plot(xValues, yValues, 'r-', 'LineWidth', 1.5); % Plot the trend line within the limited range
    
    line([9 9], ylim, 'Color', 'r', 'LineStyle', '--'); % Red vertical line at x = 8
    
    hold off;
    xlabel('Quadrants');
    ylabel('Frequency');
    title('Distribution of Waves in Hit Trials');
    grid on;
    
    % Adjust layout for better appearance
    sgtitle('Wave Analysis for Hit Trials');  % Overall title for the figure
end
