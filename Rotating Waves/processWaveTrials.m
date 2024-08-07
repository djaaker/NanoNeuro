function waves = processWaveTrials(cueTrace, parameters)

    % if strcmp(hitOrMiss, 'miss')
    %     cueTrace = 'cueMissTrace';
    % elseif strcmp(hitOrMiss, 'hit')
    %     cueTrace = 'cueTrace';
    % end

    nTrials = numel(cueTrace); 
%     waves = repmat(struct('waveTime', [], 'wavePresent', [], 'source', [], ...
%                     'rho', [], 'angVel', [], 'curl', [], 'ang', [], ...
%                     'start', [], 'end', [], 'duration', [], 'quadrant', [], 'cue', []),1,nTrials);

    waves = struct();
    
    if nTrials > 1
        nTrials = 3;
        for trial = 1:nTrials
%         for trial = 1:3
                           
            %identities
            xf = cueTrace(trial).xf;
            % trialWaves = struct();
                
            [rows,cols,~] = size( xf ); 
            [X,Y] = meshgrid( 1:cols, 1:rows ); % create grid
            
            raw = cueTrace(trial).rawLFP;
            
            % wideband filter
            xf = bandpass_filter( raw, parameters.f(1), parameters.f(2), ...
                parameters.filter_order, parameters.Fs );
            
            % GP representation
            xgp = generalized_phase( xf, parameters.Fs, parameters.lp );
            % p = xgp(:,:,start_time:stop_time); 
            p = xgp;    
            evaluation_points = find_evaluation_points( p, parameters.evaluation_angle, parameters.tol );
            
            % calculate phase gradient
            [pm,pd,dx,dy] = phase_gradient_complex_multiplication( p, parameters.pixel_spacing );
            
            %% For Rotating Waves
            % compute the source point of any detected rotating waves
            [source_rot, cav, curlz] = computeRotationalSouces(p,dx,dy,evaluation_points,X,Y);
            
            % phase correlation across polar planes
            rho_rot = NaN( 1, length(evaluation_points) ); 
        
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
                if rho_rot(jj) > parameters.plot_rho_value  
                    % find possible waves using the 40ms around the source and evaluation points
                    [wavePossibleRho(:,jj), wavePossible(:,jj), sourceCurrent(jj,1), sourceCurrent(jj,2)] = findPossibleWaves(source_rot, evaluation_points, p, curlz, jj);
                end
            end
        
            sourceCurrent(any(isnan(sourceCurrent), 2), :) = [];
            wavePossible(:, all(isnan(wavePossible))) = []; 
            wavePossibleRho(:, all(isnan(wavePossibleRho))) = [];
        
        
            [waveTime, waveRho] = extractWave(wavePossible, wavePossibleRho, parameters.plot_rho_value);
            waveTime = waveTime(~cellfun(@isempty, waveTime));
            nWaves = size(waveTime,2);
        
            % initialize all 1D arrays
            angVel = 1:nWaves; ang = angVel; waveCurl = angVel; 
    %         startTime = angVel; endTime = angVel; 
    %         duration = angVel; quadrant = angVel; cue = angVel; iFreq = angVel; 
        
            for jj = 1:nWaves
                if ~isempty(waveTime{jj})
                    % characterize the rotation of the wave using curl and angular velocity
                    [angVel(jj), ang(jj), waveCurl(jj)] = rotationalCharacteristics(waveTime{jj}, cav, curlz, sourceCurrent(jj,:));
                end
            end
        
            allTimes = cell2mat(cellfun(@(x) x(:)', waveTime(~cellfun(@isempty, waveTime)), 'UniformOutput', false));
            wavePresent = zeros(1,3001);
            wavePresent(allTimes) = 1;
        
            waves(trial).waveTime = waveTime;
            waves(trial).wavePresent = wavePresent;
            waves(trial).source =  sourceCurrent;
            waves(trial).rho = waveRho;
            waves(trial).angVel = angVel;
            waves(trial).curl =  waveCurl;
            waves(trial).ang = ang;
            
            waves(trial).start = cellfun(@(x) x(1), waveTime(~cellfun(@isempty, waveTime)));
            waves(trial).end = cellfun(@(x) x(end), waveTime(~cellfun(@isempty, waveTime)));
            waves(trial).duration = waves(trial).end - waves(trial).start;
            
            waves(trial).quadrant = ceil((waves(trial).start / 3001) * 16);
            waves(trial).cue = floor(waves(trial).start / 1500);
                   
            if parameters.verbose == true
    %             fprintf('trial: %d, ',trial)
                if mod(trial,10) == 0
                    fprintf("trial %d done\n",trial)
                end
            end
%             fprintf('trial %d\nwaves(trial).start %d\n', trial, waves(trial).start)
        end % end of trial analysis 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        else
            %identities
            xf = cueTrace.xf;
            % trialWaves = struct();
                
            [rows,cols,~] = size( xf ); 
            [X,Y] = meshgrid( 1:cols, 1:rows ); % create grid
            
            raw = cueTrace.rawLFP;
            
            % wideband filter
            xf = bandpass_filter( raw, parameters.f(1), parameters.f(2), ...
                parameters.filter_order, parameters.Fs );
    
            disp(xf(4,4,1))
            
            % GP representation
            xgp = generalized_phase( xf, parameters.Fs, parameters.lp );
            % p = xgp(:,:,start_time:stop_time); 
            p = xgp;    
            evaluation_points = find_evaluation_points( p, parameters.evaluation_angle, parameters.tol );
            
            % calculate phase gradient
            [pm,pd,dx,dy] = phase_gradient_complex_multiplication( p, parameters.pixel_spacing );
            
            %% For Rotating Waves
            % compute the source point of any detected rotating waves
            [source_rot, cav, curlz] = computeRotationalSouces(p,dx,dy,evaluation_points,X,Y);
            
            % phase correlation across polar planes
            rho_rot = NaN( 1, length(evaluation_points) ); 
        
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
                if rho_rot(jj) > parameters.plot_rho_value  
                    % find possible waves using the 40ms around the source and evaluation points
                    [wavePossibleRho(:,jj), wavePossible(:,jj), sourceCurrent(jj,1), sourceCurrent(jj,2)] = findPossibleWaves(source_rot, evaluation_points, p, curlz, jj);
                end
            end
        
            sourceCurrent(any(isnan(sourceCurrent), 2), :) = [];
            wavePossible(:, all(isnan(wavePossible))) = []; 
            wavePossibleRho(:, all(isnan(wavePossibleRho))) = [];
        
        
            [waveTime, waveRho] = extractWave(wavePossible, wavePossibleRho, parameters.plot_rho_value);
            waveTime = waveTime(~cellfun(@isempty, waveTime));
            nWaves = size(waveTime,2);
        
            % initialize all 1D arrays
            angVel = 1:nWaves; ang = angVel; waveCurl = angVel; 
    %         startTime = angVel; endTime = angVel; 
    %         duration = angVel; quadrant = angVel; cue = angVel; iFreq = angVel; 
        
            for jj = 1:nWaves
                if ~isempty(waveTime{jj})
                    % characterize the rotation of the wave using curl and angular velocity
                    [angVel(jj), ang(jj), waveCurl(jj)] = rotationalCharacteristics(waveTime{jj}, cav, curlz, sourceCurrent(jj,:));
                end
            end
        
            allTimes = cell2mat(cellfun(@(x) x(:)', waveTime(~cellfun(@isempty, waveTime)), 'UniformOutput', false));
            wavePresent = zeros(1,3001);
            wavePresent(allTimes) = 1;
        
            waves.waveTime = waveTime;
            waves.wavePresent = wavePresent;
            waves.source =  sourceCurrent;
            waves.rho = waveRho;
            waves.angVel = angVel;
            waves.curl =  waveCurl;
            waves.ang = ang;
            
            waves.start = cellfun(@(x) x(1), waveTime(~cellfun(@isempty, waveTime)));
            waves.end = cellfun(@(x) x(end), waveTime(~cellfun(@isempty, waveTime)));
            waves.duration = waves.end - waves.start;
            
            waves.quadrant = ceil((waves.start / 3001) * 16);
            waves.cue = floor(waves.start / 1500);
%             fprintf('trial %d\nangVel %d\nwaves.start %d\n', trial, angVel, waves.start)
    end
end     

