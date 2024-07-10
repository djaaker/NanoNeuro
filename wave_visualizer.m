%function wave_visualizer(behaviour,trial)

% parameters
parameters.Fs = 1000; % data sampling rate [Hz]
parameters.filter_order = 4; parameters.f = [5 20]; % filter parameters
parameters.start_time = -300; % time point to start analysis
parameters.stop_time = 0; % stop analysis (relative to target onset)
parameters.pixel_spacing = 0.4; % spacing between electrodes [mm]
parameters.lp = 0; % cutoff for negative frequency detection [Hz]
parameters.evaluation_angle = pi; parameters.tol = 0.2; % evaluation points
parameters.windowBeforeCue = 1.5; % in seconds 
parameters.windowAfterCue = 1.5; % in seconds 
plot_rho_value = 0.7; plot_time = 20; pause_length = 0.2; 
parameters.extract = true;


 
if parameters.extract == true 
    rotWaves = struct();
    cnt = 1;
    wavesHit = repmat(struct('waveTime', [], 'wavePresent', [], 'source', [], ...
                    'rho', [], 'angVel', [], 'curl', [], 'ang', [], ...
                    'start', [], 'end', [], 'duration', [], 'quadrant', [], 'cue', []),1,143);

for trial = 1:numel(behaviour.cueHitTrace)
% for trial = 1:3
% trial = 1;
    
    
    
    %identities
    xf = behaviour.cueHitTrace(trial).xf;
    % trialWaves = struct();
    
    % Initialize trialWaves structure with empty arrays or appropriate initial values
    
    [rows,cols,~] = size( xf ); 
    [X,Y] = meshgrid( 1:cols, 1:rows ); % create grid
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    raw = behaviour.cueHitTrace(trial).rawLFP;
    
    % wideband filter
    xf = bandpass_filter( raw, parameters.f(1), parameters.f(2), ...
        parameters.filter_order, parameters.Fs );
    
    % GP representation
    xgp = generalized_phase( xf, parameters.Fs, parameters.lp );
    % p = xgp(:,:,start_time:stop_time); 
    p = xgp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %p = behaviour.cueHitTrace(trial).xgp;
    
    evaluation_points = find_evaluation_points( p, parameters.evaluation_angle, parameters.tol );
    
    % calculate phase gradient
    [pm,pd,dx,dy] = phase_gradient_complex_multiplication( p, parameters.pixel_spacing );
    
    %% For Expanding Waves
    % divergence calculation
    source = find_source_points( evaluation_points, X, Y, dx, dy );
    
    % % phase correlation with distance (\rho_{\phi,d} measure)
    % rho_exp = zeros( 1, length(evaluation_points));
    % for jj = 1:length(evaluation_points)
    %     ph = angle( p(:,:,evaluation_points(jj)) );
    %     rho_exp(jj) = phase_correlation_distance( ph, source(:,jj), parameters.pixel_spacing );
    % end

    %% For Rotating Waves
    % compute the source point of any detected rotating waves
    [source_rot, cav, curlz] = computeRotationalSouces(p,dx,dy,evaluation_points,X,Y);
    
    % phase correlation across polar planes
    rho_rot = NaN( 1, length(evaluation_points) ); wave_cnt = 1;

    % initialize 2D arrays
    sourceCurrent = NaN(length(evaluation_points),2); 

    % initialize variable arrays
    wavePossibleRho = NaN(41,length(evaluation_points)); wavePossible = wavePossibleRho;    

    for jj = 1:length(evaluation_points)
        
        pl = angle( p(:,:,evaluation_points(jj)) );
        rho_rot(jj) = abs(phase_correlation_rotation( pl, curlz(:,:,jj), source_rot(:,jj)));

    end
    
    for jj = 1:length(rho_rot)
        % if there is a rotating wave, capture just the wave
        if rho_rot(jj) > plot_rho_value  
            % find possible waves using the 40ms around the source and evaluation points
            [wavePossibleRho(:,jj), wavePossible(:,jj), sourceCurrent(jj,1), sourceCurrent(jj,2)] = findPossibleWaves(source_rot, evaluation_points, p, curlz, jj);
        end
    end

    sourceCurrent(any(isnan(sourceCurrent), 2), :) = [];
    wavePossible(:, all(isnan(wavePossible))) = []; 
    wavePossibleRho(:, all(isnan(wavePossibleRho))) = [];


    [waveTime, waveRho] = extractWave(wavePossible, wavePossibleRho, plot_rho_value);
    waveTime = waveTime(~cellfun(@isempty, waveTime));
    nWaves = size(waveTime,2);

    % initialize all 1D arrays
    angVel = 1:nWaves; ang = angVel; waveCurl = angVel; startTime = angVel; endTime = angVel; 
    duration = angVel; quadrant = angVel; cue = angVel; iFreq = angVel; 

    for jj = 1:nWaves
        if ~isempty(waveTime{jj})
            % characterize the rotation of the wave using curl and angular velocity
            [angVel(jj), ang(jj), waveCurl(jj)] = rotationalCharacteristics(waveTime{jj}, cav, curlz, sourceCurrent(jj,:));
        end
    end

    allTimes = cell2mat(cellfun(@(x) x(:)', waveTime(~cellfun(@isempty, waveTime)), 'UniformOutput', false));
    wavePresent = zeros(1,3001);
    wavePresent(allTimes) = 1;

    wavesHit(cnt).waveTime = waveTime;
    wavesHit(cnt).wavePresent = wavePresent;
    wavesHit(cnt).source =  sourceCurrent;
    wavesHit(cnt).rho = waveRho;
    wavesHit(cnt).angVel = angVel;
    wavesHit(cnt).curl =  waveCurl;
    wavesHit(cnt).ang = ang;
    
    wavesHit(cnt).start = cellfun(@(x) x(1), waveTime(~cellfun(@isempty, waveTime)));
    wavesHit(cnt).end = cellfun(@(x) x(end), waveTime(~cellfun(@isempty, waveTime)));
    wavesHit(cnt).duration = wavesHit(cnt).end - wavesHit(cnt).start;
    
    wavesHit(cnt).quadrant = ceil((wavesHit(cnt).start / 3001) * 16);
    wavesHit(cnt).cue = floor(wavesHit(cnt).start / 1500);
    
    cnt = cnt + 1;
    % wavesHit(trial) = trialWaves;

    % disp(wavesHit(cnt).ang(1))

    fprintf('trial: %d, ',trial)
    if mod(trial,10) == 0
        fprintf("trial %d done\n",trial)
    end
end
end     % end of trial analysis 

rotWaves.wavesHit = wavesHit;

% end  

% save('RotWavesHits_40-80.mat',"rotWaves")

% rotWaves = load('RotWavesHits_5-40');
hitsVel = []; hitsStart = []; hitsQuad = []; wavesPresAll = NaN(3001,143);
for i = 1:numel(behaviour.cueHitTrace)
    % wavesPresAll(:,i) = rotWaves.wavesHit(i).wavePresent;
    hitsVel = [hitsVel rotWaves.wavesHit(i).angVel];
    hitsStart = [hitsStart rotWaves.wavesHit(i).start];
    hitsQuad = [hitsQuad rotWaves.wavesHit(i).quadrant];
end

plotWaveRaster(rotWaves.wavesHit,[],behaviour,[])



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
hold off;

% Add a red vertical line in the center
xCenter = (max(hitsStart) + min(hitsStart)) / 2;
line([xCenter xCenter], get(gca, 'ylim'), 'Color', 'r', 'LineWidth', 1.5);

% Distribution of waves in hit trials
subplot(2, 1, 2);
histogram(hitsQuad, 'BinLimits', [0 16], 'BinWidth', 1);
hold on;
xlim([0 16]);

% Fit a normal distribution to the data
pd = fitdist(hitsQuad', 'Normal');
xValues = linspace(0, 16, 100); % X values from 0 to 16
yValues = pdf(pd, xValues) * numel(hitsQuad) * (xValues(2) - xValues(1)); % Scaled PDF

plot(xValues, yValues, 'r-', 'LineWidth', 1.5); % Plot the trend line within the limited range

line([8 8], ylim, 'Color', 'r', 'LineStyle', '--'); % Red vertical line at x = 8

hold off;
xlabel('Quadrants');
ylabel('Frequency');
title('Distribution of Waves in Hit Trials');
grid on;

% Adjust layout for better appearance
sgtitle('Wave Analysis for Hit Trials');  % Overall title for the figure

