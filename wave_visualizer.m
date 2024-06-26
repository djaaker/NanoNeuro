%function pca_wave(behaviour,trial)

%identities
behaviour = IntanBehaviourBaseline;

xf = behaviour.cueHitTrace(trial).xf;

% options
options.subject = 'W'; % this can be 'W' or 'T' (two marmoset subjects)
options.plot = true; % this option turns plots ON or OFF
options.plot_shuffled_examples = false; % example plots w/channels shuffled in space

% parameters
parameters.Fs = 1000; % data sampling rate [Hz]
parameters.filter_order = 4; parameters.f = [5 40]; % filter parameters
parameters.start_time = -300; % time point to start analysis
parameters.stop_time = 0; % stop analysis (relative to target onset)
parameters.pixel_spacing = 0.4; % spacing between electrodes [mm]
parameters.lp = 0; % cutoff for negative frequency detection [Hz]
parameters.evaluation_angle = pi; parameters.tol = 0.2; % evaluation points
plot_rho_value = 0.5; plot_time = 20; pause_length = 0.2; 


[rows,cols,~] = size( xf ); 
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

ctr = 1;
for jj = 1:length(evaluation_points)
    
    % animated wave plot
    if ( rho(jj) > plot_rho_value )
        
        % get start and stop time, truncating if necessary
        st = evaluation_points(jj) - plot_time; sp = evaluation_points(jj) + plot_time;
        if ( st < 1 ), st = 1; end; if ( sp > size(xf,3) ), sp = size(xf,3); end

        % get data to plot, shuffle if option is chosen
        x_plot = xf(:,:,st:sp); 
    end
end
