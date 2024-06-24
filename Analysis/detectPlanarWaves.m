function [Waves] = detectPlanarWaves(xf,xgp,wt,behaviourTrace,parameters,threshold)

spacing = parameters.spacing;
pgdThres = threshold;
dirThres = 15; %in degrees
X = parameters.X;
Y = parameters.Y;

waveTimeWindow = 40; % in points

for ii=1:size(behaviourTrace,2)
    Waves(ii).xf = xf{1,ii};
    Waves(ii).trialTime = behaviourTrace(ii).LFPIndex;
    Waves(ii).p = xgp{1,ii};
    Waves(ii).wt = wt{1,ii};
    %p = arrayfun(@(jj) inpaint_nans(p(:,:,jj)),1:size(p,3));
    Waves(ii).evaluationPoints = find_evaluation_points(Waves(ii).p,pi,0.2);
%     plot_evaluation_points( Waves(ii).p, Waves(ii).evaluationPoints );
    [Waves(ii).pm,Waves(ii).pd,Waves(ii).dx,Waves(ii).dy] = phase_gradient_complex_multiplication( Waves(ii).p, spacing );
    % Phase gradient directionality 
    [Waves(ii).PGD] = phase_gradient_directionality(Waves(ii).pm,Waves(ii).dx,Waves(ii).dy);
    % Wavelength
    Waves(ii).wl = 1./abs(Waves(ii).pm);
    % Instantaneous speed
    Waves(ii).insts = instantaneous_speed(Waves(ii).wt,Waves(ii).pm);
    Waves(ii).s = speedSpatial(Waves(ii).wt,Waves(ii).pm);
    Waves(ii).l = wavelengthSpatial(Waves(ii).pm);
    % Velcity direction unit vector
    [Waves(ii).vx, Waves(ii).vy] = wavefront_direction(Waves(ii).pd,Waves(ii).insts);
    Waves(ii).velDir = atan2(Waves(ii).vy,Waves(ii).vx);
    % divergence calculation
    Waves(ii).source = find_source_points( Waves(ii).evaluationPoints, X, Y, Waves(ii).dx, Waves(ii).dy );
    % phase correlation with distance (\rho_{\phi,d} measure)
    Waves(ii).rho = zeros( 1, length(Waves(ii).evaluationPoints) );
    for jj = 1:length(Waves(ii).evaluationPoints)
        st = Waves(ii).evaluationPoints(jj)-waveTimeWindow;
        if st<=0, st = 1; end  % adjusting is waves is detected at the start of window
        sp = Waves(ii).evaluationPoints(jj)+waveTimeWindow;
        if sp>size(Waves(ii).p,3), sp = Waves(ii).evaluationPoints(jj); end % adjusting is waves is detected at the end of window
        pgd = Waves(ii).PGD(1,st:sp);
        dir = Waves(ii).velDir(1,st:sp);
        waveClusters = findThresSeg(pgd,pgdThres);
        % Rejecting segements less than 9ms
        clusterSize = diff(waveClusters,1,2);
        rejectIndex = find(clusterSize<9);
        waveClusters(rejectIndex,:) = [];
        % Getting segments with cummulative sum of first differences less
        % than 15 degress
        for aa=1:size(waveClusters,1)
            dirSegments = segmentDirArray(dir(waveClusters(aa,1):waveClusters(aa,2)),15);
            dirSegSize = diff(dirSegments,1,2);
            [maxvalue, maxindex] = max(dirSegSize);
            %Rejecting segments that that less than 9 seconds long
            if maxvalue<=9
                waveClusters(aa,:) = NaN;
                continue  
            end
            %Updating the waveCluster start and stop points
            waveClusters(aa,1) = waveClusters(aa,1)+dirSegments(maxindex,1)-1;
            waveClusters(aa,2) = waveClusters(aa,1)+dirSegments(maxindex,2)-dirSegments(maxindex,1);
        end
        waveClusters = removeNaNRows(waveClusters);
        % Selecting the wavecluster with maximum duration 
        % This is done to have just one wave at each evaluation point
        if isempty(waveClusters)
            Waves(ii).evaluationPoints(jj) = NaN;
            Waves(ii).waveTime(jj,:) = [NaN,NaN];
            continue
        end
        clusterSize = diff(waveClusters,1,2);
        [~, clusterMaxIndx] = max(clusterSize);
        Waves(ii).waveTime(jj,:) = [st-1+waveClusters(clusterMaxIndx,1),st-1+waveClusters(clusterMaxIndx,2)];
        Waves(ii).evaluationPoints(jj) = findArrayMidPoint([Waves(ii).waveTime(jj,1):Waves(ii).waveTime(jj,2)]);
    end
    % Removing evaluation points in which no wave was detected 
    Waves(ii).evaluationPoints = Waves(ii).evaluationPoints(~isnan(Waves(ii).evaluationPoints));
    Waves(ii).waveTime = removeNaNRows(Waves(ii).waveTime);
    % Removing duplicate detected waves
%     Waves(ii).evaluationPoints = unique((Waves(ii).evaluationPoints)','rows','stable');
%     Waves(ii).waveTime = unique(Waves(ii).waveTime,'rows','stable');
    % Merging repetetive waves 
    Waves(ii).nWaves = size(Waves(ii).waveTime,1);
    kk = 1;
    while kk < Waves(ii).nWaves
        if(Waves(ii).waveTime(kk,2) > Waves(ii).waveTime(kk+1,1)) % check there is a overlap
            if (Waves(ii).waveTime(kk,2) < Waves(ii).waveTime(kk+1,2))
                Waves(ii).waveTime(kk,2) = Waves(ii).waveTime(kk+1,2); % Merge the two waves
            end
            Waves(ii).evaluationPoints(kk) = Waves(ii).waveTime(kk,1);
            Waves(ii).waveTime(kk+1,:) = [];
            Waves(ii).evaluationPoints(kk+1) = [];
            Waves(ii).nWaves = Waves(ii).nWaves-1;
        else
            kk = kk+1;
        end
    end
    
%     plot_evaluation_points( Waves(ii).p, Waves(ii).evaluationPoints );
    % Calculating the speed, amplitude and wave direction of the detected waves in
    % that trial window
    for kk = 1:Waves(ii).nWaves
        %Waves(ii).speed(kk) = mean(abs(Waves(ii).insts(:,:,st:sp)),[1 2 3]); % speed in cm/s
        Waves(ii).speed(kk) = mean(abs(Waves(ii).s(Waves(ii).waveTime(kk,1):Waves(ii).waveTime(kk,2))),'all'); % speed in cm/s
        Waves(ii).waveDir(kk) = mean(Waves(ii).velDir(Waves(ii).waveTime(kk,1):Waves(ii).waveTime(kk,2)),'all'); 
        Waves(ii).wavelength(kk) = mean(Waves(ii).l(Waves(ii).waveTime(kk,1):Waves(ii).waveTime(kk,2)),'all');
        Waves(ii).waveDuration(kk) = Waves(ii).waveTime(kk,2)-Waves(ii).waveTime(kk,1)+1;
        Waves(ii).waveAmp(kk) = max(abs(Waves(ii).p(:,:,Waves(ii).waveTime(kk,1):Waves(ii).waveTime(kk,2))),[],"all");
    end
    %Waves(ii).speedpdg = pgdMean(Waves(ii).PGD,Waves(ii).s,0.51);
    %Waves(ii).dirpdg = pgdMean(Waves(ii).PGD,Waves(ii).velDir,0.51);
    %plot_vector_field( exp( 1i .* Waves(ii).pd(:,:,100) ), 0 );
    %plot_wave_examples( LFP.xf(:,:,Encoder.trialTime(ii,3):Encoder.trialTime(ii,4)), options, ii, Waves(ii).evaluationPoints, Waves(ii).source, Waves(ii).rho );
end
  
end

