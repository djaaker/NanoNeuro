function [wavePossibleRho, wavePossible, xSource, ySource] = findPossibleWaves(source_rot, evaluation_points,p, curlz, jj)

sourceCurrent = source_rot(:, jj);

% cycle through each time point before/after the eval point to
% determine if it's a wave
for kk = -20:20
    mid_point = evaluation_points(jj);
    % fprintf("\neval_points(jj) = %d\nkk + 21 = %d\nstart + kk = %d\n",evaluation_points(jj),kk+21,start+kk)
    
    % start point for wave detected in first 20ms is set to the start of trial
    if evaluation_points(jj) + kk < 1
        continue; 
    else

        % fprintf("\nmidpoint: %d\nkk: %d\n", mid_point, kk) 

        wp_st = mid_point - 20; wp_sp = mid_point + 20;

        if wp_st < 1
            wp_st = 1;
        end
        if wp_sp > 3001
            wp_sp = 3001;
        end
        
        wavePossible = wp_st : wp_sp; t = mid_point + kk + 21;
        pl_poss = nan(8, 8, numel(wavePossible)); 
        
        if t > 3001
            continue;
        end

        wavePossibleRho = nan(size(wavePossible)); 
        for qq = 1:length(wavePossible)
            pl_poss(:, :, qq) = angle(p(:, :, qq));
            wavePossibleRho(qq) = abs(phase_correlation_rotation(pl_poss(:, :, qq), curlz(:, :, t), sourceCurrent));
        end

        xSource = sourceCurrent(1); ySource = sourceCurrent(2);       
    end
end

if length(wavePossible) < 41
    numZeros = 41 - length(wavePossible);
    wavePossibleRho = [zeros(1, numZeros), wavePossibleRho]; 
    wavePossible = [zeros(1, numZeros), wavePossible];
end            
    
end
