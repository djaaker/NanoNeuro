function wave_visualizer(behaviour,waves,trial,wave_id)


% options
options.plot = true; % this option turns plots ON or OFF
options.plot_shuffled_examples = false; % example plots w/channels shuffled in space
options.subject = 'W';

% parameters
parameters.Fs = 1000; % data sampling rate [Hz]
parameters.filter_order = 4; parameters.f = [5 40]; % filter parameters
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
start_ids = find(waves.wavesHit(trial).waveStart==1); start_time = start_ids(wave_id);
duration = waves.wavesHit(trial).waveDuration(wave_id);
wave_points =  start_time : start_time + duration;


[rows,cols,~] = size( behaviour.cueHitTrace(trial).xf ); channels = rows*cols;
[X,Y] = meshgrid( 1:cols, 1:rows ); % create grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = behaviour.cueHitTrace(trial).xgp;

evaluation_points = find_evaluation_points( p, parameters.evaluation_angle, parameters.tol );

% calculate phase gradient
[pm,pd,dx,dy] = phase_gradient_complex_multiplication( p, parameters.pixel_spacing );

% divergence calculation
source = find_source_points( evaluation_points, X, Y, dx, dy );

% phase correlation with distance (\rho_{\phi,d} measure)
rho = zeros( 1, length(evaluation_points) );
for jj = 1:length(evaluation_points)
    
    ph = angle( p(:,:,evaluation_points(jj)) );
    rho(jj) = phase_correlation_distance( ph, source(:,jj), parameters.pixel_spacing );
    
end

% plotting 2 - wave examples
if options.plot, plot_wave_examples( xf(:,:,wave_points), options, trial, evaluation_points, source, rho ); end

for jj = 1:length(evaluation_points)
    
    % animated wave plot
    if ( rho(jj) > plot_rho_value )
        
        % get start and stop time, truncating if necessary
        st = evaluation_points(jj) - plot_time; sp = evaluation_points(jj) + plot_time;
        if ( st < 1 ), st = 1; end; if ( sp > size(xf,3) ), sp = size(xf,3); end

        % get data to plot, shuffle if option is chosen
        x_plot = xf(:,:,st:sp); 
        if ( options.plot_shuffled_examples == true ), x_plot = shuffle_channels( x_plot ); end        
        
        % create plot
        figure; title( sprintf( 'trial %d, wave example %d, 0 of %d ms', trial, ctr, size(x_plot,3) ) );
        color_range = [ min(reshape(x_plot,[],1)) max(reshape(x_plot,[],1)) ];
        h = imagesc( x_plot(:,:,1) ); hold on; axis image; 
        plot( source(1,jj), source(2,jj), '.', 'markersize', 35, 'color', [.7 .7 .7] );
        set( gca, 'linewidth', 3, 'xtick', [], 'ytick', [], 'fontname', 'arial', 'fontsize', 16, 'ydir', 'reverse' ); 
        colormap( M.myMap ); box on; xlabel( 'electrodes' ); ylabel( 'electrodes' ); clim( color_range )
        
        % create colorbar
        cb = colorbar();
    
        set( cb, 'location', 'southoutside' )
        set( cb, 'position', [0.6661    0.1674    0.2429    0.0588] );

        set( get(cb,'ylabel'), 'string', 'Amplitude (\muV)' ); set( cb, 'linewidth', 2 )
        
        % animate plot
        for kk = 1:size(x_plot,3)
            set( h, 'cdata', x_plot(:,:,kk) ); 
            set( get(gca,'title'), 'string', ...
                sprintf( 'trial %d, wave example %d, %d of %d ms', trial, ctr, kk, size(x_plot,3) ) )
            pause(pause_length); 
        end
    end
end

v = VideoWriter('wave_1.avi');

% plotting 1 - evaluation points
%if options.plot, plot_evaluation_points( p, evaluation_points ); end


end