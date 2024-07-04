%function wave_visualizer(behaviour,trial)


% options
options.plot = false; % this option turns plots ON or OFF
options.plot_shuffled_examples = false; % example plots w/channels shuffled in space
options.subject = 'W';
trial = 1;

% init
M = load( 'generalized-phase-main\input\myMap.mat' );

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
plot_rho_value = 0.5; plot_time = 20; pause_length = 0.2; 

%identities
xf = behaviour.cueHitTrace(trial).xf;
%start_ids = find(waves.wavesHit(trial).waveStart==1); 
%duration = waves.wavesHit(trial).waveDuration(wave_id);
%wave_points =  start_time : start_time + duration;


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

% phase correlation with distance (\rho_{\phi,d} measure)
rho_exp = zeros( 1, length(evaluation_points) );
for jj = 1:length(evaluation_points)
    
    ph = angle( p(:,:,evaluation_points(jj)) );
    rho_exp(jj) = phase_correlation_distance( ph, source(:,jj), parameters.pixel_spacing );
    
end


%% For Rotating Waves
rotWaves = {}; cav = nan(size(p));

% curl calculation
for ii = 1:length(xf(1,1,:))
    dx(:,:,ii) = inpaint_nans(dx(:,:,ii));
    dy(:,:,ii) = inpaint_nans(dy(:,:,ii));
    
    %cav(:,:,ii) = curl(X,Y,dx(:,:,evaluation_points(ii)), dy(:,:,evaluation_points(ii)));
    cav(:,:,ii) = curl(X,Y,dx(:,:,ii), dy(:,:,ii));
    cav(:,:,ii) = inpaint_nans(cav(:,:,ii));
end


% source point
source_rot = nan( 2, length(evaluation_points) );
for ii = 1:length(evaluation_points)    
    [yy,xx] = find( cav(:,:,ii) == max( reshape(cav(:,:,ii), 1, [] ) ) );
    if numel(yy) == 1
        source_rot(1,ii) = xx; source_rot(2,ii) = yy;
    else
        source_rot(1,ii) = NaN; source_rot(2,ii) = NaN; 
    end
end


% phase correlation across polar planes
rho_rot = zeros( 1, length(evaluation_points) ); wave_cnt = 1;
for jj = 1:length(evaluation_points)
    
    pl = angle( p(:,:,evaluation_points(jj)) ); waveTime = [];
    rho_rot(jj) = abs(phase_correlation_rotation( pl, cav(:,:,jj), source_rot(:,jj)));

    % if there is a rotating wave, capture just the wave
    if rho_rot(jj) > plot_rho_value
        % preallocate structures and define source point
         
        source_current = source_rot(:, jj);
        
        % cycle through each time point before/after the eval point to
        % determine if it's a wave
        for kk = -20:20
            mid_point = evaluation_points(jj);
            % fprintf("\neval_points(jj) = %d\nkk + 21 = %d\nstart + kk = %d\n",evaluation_points(jj),kk+21,start+kk)
            
            % start point for wave detected in first 20ms is set to the start of trial
            if evaluation_points(jj) + kk < 1
                continue; 
            else

                fprintf("\nmidpoint: %d\nkk: %d\n", mid_point, kk) 

                wp_st = mid_point - 20; wp_sp = mid_point + 20;
    
                if wp_st < 1
                    wp_st = 1;
                end
                if wp_sp > 3001
                    wp_sp = 3001;
                end
                
                wave_possible = wp_st : wp_sp; t = mid_point + kk + 21;
                pl_poss = nan(8, 8, numel(wave_possible)); 
                
                if t > 3001
                    continue;
                end
                
                wave_possible_rho = nan(size(wave_possible)); 
                for qq = 1:length(wave_possible)
                    pl_poss(:, :, qq) = angle(p(:, :, qq));
                    wave_possible_rho(qq) = abs(phase_correlation_rotation(pl_poss(:, :, qq), cav(:, :, t), source_current));
                end
                
                % wave_possible_rho = wave_possible_rho(rho_real);
            
            end
        end
        
        % find the time that the wave is present
        waveTime = extract_wave(wave_possible, wave_possible_rho);
        fprintf("\nwave number:%d\nwave midpoint: %d\nwaveTime start: %d\n",wave_cnt, mid_point, waveTime(1))
    
        % populate the rotWaves struct
        rotWaves(wave_cnt).waveTime = waveTime; 
        rotWaves(wave_cnt).source = source_current;
        rotWaves(wave_cnt).wavePossible = wave_possible;
        rotWaves(wave_cnt).rho = wave_possible_rho;
        % rotWaves.angVel = cav(:, :, wave_possible);
    
        wave_cnt = wave_cnt + 1;
    end

end


eval_pts_rot = []; eval_pts_exp = [];
for ii = 1:length(evaluation_points)
    if rho_rot(ii) > rho_exp(ii)
        eval_pts_rot = [eval_pts_rot evaluation_points(ii)];
    else
        eval_pts_exp = [eval_pts_exp evaluation_points(ii)];
    end

end 


% Initialize video writer object
% videoFileName = 'rot_waves_80_200.avi'; 
% v = VideoWriter(videoFileName, 'Motion JPEG AVI'); v.FrameRate = 5;
% open(v);

if options.plot == true
    ctr = 1;
    for jj = 1:numel(rotWaves)
        % animate plot
        tw = rotWaves(jj); time = tw.waveTime;
        if ~isempty(time) && (max(tw.rho) > .75)
            st = time(1); sp = time(end);
            wave_plot = xf(:,:,st:sp);
            
            %create plot
            figure; title( sprintf( 'trial %d, wave example %d, 0 of %d ms', trial, ctr, size(wave_plot,3) ) );
            color_range = [ min(reshape(wave_plot,[],1)) max(reshape(wave_plot,[],1)) ]; 
            h = imagesc( wave_plot(:,:,1) ); hold on; axis image; 
            plot( source(1,jj), source(2,jj), '.', 'markersize', 35, 'color', [.7 .7 .7] );
            set( gca, 'linewidth', 3, 'xtick', [], 'ytick', [], 'fontname', 'arial', 'fontsize', 16, 'ydir', 'reverse' ); 
            colormap( M.myMap ); box on; xlabel( 'electrodes' ); ylabel( 'electrodes' ); clim( color_range )
    
            % create colorbar
            cb = colorbar();
    
            set( cb, 'location', 'southoutside' )
            set( cb, 'position', [0.6661    0.1674    0.2429    0.0588] );
    
            set( get(cb,'ylabel'), 'string', 'Amplitude (\muV)' ); set( cb, 'linewidth', 2 )
    
            % animate plot
            for kk = 1:size(wave_plot,3)
                set( h, 'cdata', wave_plot(:,:,kk) ); 
                set( get(gca,'title'), 'string', ...
                    sprintf( 'trial %d, wave example %d, %d of %d ms', trial, ctr, kk, size(wave_plot,3) ) )
                pause(pause_length); 
                % frame = getframe(gcf);
                % writeVideo(v,frame)
            end
    
            ctr = ctr + 1;
            %close(gcf);

        end
    end
end

    
        % % animated wave plot
        % %if ( rho_rot(jj) > plot_rho_value )
        % 
        %     % get start and stop time, truncating if necessary
        %     st = eval_pts_rot(jj) - plot_time; sp = eval_pts_rot(jj) + plot_time;
        %     if ( st < 1 ), st = 1; end; if ( sp > size(xf,3) ), sp = size(xf,3); end
        % 
        %     % get data to plot, shuffle if option is chosen
        %     x_plot = xf(:,:,st:sp); 
        %     if ( options.plot_shuffled_examples == true ), x_plot = shuffle_channels( x_plot ); end        
        % 
        %     % create plot
        %     figure; title( sprintf( 'trial %d, wave example %d, 0 of %d ms', trial, ctr, size(x_plot,3) ) );
        %     color_range = [ min(reshape(x_plot,[],1)) max(reshape(x_plot,[],1)) ]; 
        %     h = imagesc( x_plot(:,:,1) ); hold on; axis image; 
        %     plot( source(1,jj), source(2,jj), '.', 'markersize', 35, 'color', [.7 .7 .7] );
        %     set( gca, 'linewidth', 3, 'xtick', [], 'ytick', [], 'fontname', 'arial', 'fontsize', 16, 'ydir', 'reverse' ); 
        %     colormap( M.myMap ); box on; xlabel( 'electrodes' ); ylabel( 'electrodes' ); clim( color_range )
        % 
        %     % create colorbar
        %     cb = colorbar();
        % 
        %     set( cb, 'location', 'southoutside' )
        %     set( cb, 'position', [0.6661    0.1674    0.2429    0.0588] );
        % 
        %     set( get(cb,'ylabel'), 'string', 'Amplitude (\muV)' ); set( cb, 'linewidth', 2 )
        % 
        %     % animate plot
        %     for kk = 1:size(x_plot,3)
        %         set( h, 'cdata', x_plot(:,:,kk) ); 
        %         set( get(gca,'title'), 'string', ...
        %             sprintf( 'trial %d, wave example %d, %d of %d ms', trial, ctr, kk, size(x_plot,3) ) )
        %         pause(pause_length); 
        %         % frame = getframe(gcf);
        %         % writeVideo(v,frame)
        %     end
        % 
        %     ctr = ctr + 1;
        %     %close(gcf);
    
       % else 
            %fprintf("\neval point %d rejected \nrho value = %0.3f < %0.3f\n", eval_pts_rot(jj), rho_rot(jj), plot_rho_value);
    
        %end
%     end
% end 
% Close the video file
%close(v);

%end