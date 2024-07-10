function [waveTimeCell, waveRhoCell] = extractWave(wavePossible, wavePossibleRho, threshold)
    if nargin < 3
        threshold = 0.5;
    end

    waveTimeCell = cell(1, size(wavePossible, 2)); 
    waveRhoCell = cell(1, size(wavePossible, 2)); 

    for jj = 1:size(wavePossible, 2)
        if ~isnan(wavePossible(1, jj))  % Check if the column of wavePossible is not NaN
            
            waveTime = [];
            waveRho = [];
            
            for i = 1:size(wavePossibleRho, 1) % Iterate over the first dimension of wavePossibleRho
                if wavePossibleRho(i, jj) > threshold
                    waveTime = [waveTime, wavePossible(i, jj)];
                    waveRho = [waveRho, wavePossibleRho(i, jj)];
                end
            end
            
            [waveTime, uniqueIndices] = unique(waveTime);
            waveRho = waveRho(uniqueIndices);
            
            if length(waveTime) > 1
                % fprintf('\nwaveTime: %d\n', waveTime(1))
                waveTime = fillGaps(waveTime);
                waveRho = updateRho(waveTime, wavePossible(:, jj), wavePossibleRho(:, jj));
                [waveTime, waveRho] = removeShortSequences(waveTime, waveRho, 9);
                if length(waveTime) > 9
                    % fprintf('    --- WAVE DETECTED ---\nTIME: %dms\n', length(waveTime))
                end
            else 
                % fprintf('\nRejected\n')
            end
            
            waveTimeCell{jj} = waveTime;  % Store the result in the cell array
            waveRhoCell{jj} = waveRho;    % Store the result in the cell array
        end
    end
end


function filledArray = fillGaps(inputArray)
    inputArray = sort(inputArray);
    filledArray = inputArray(1);

    for i = 2:length(inputArray)
        gap = inputArray(i) - inputArray(i-1);
        
        if gap > 1 && gap < 3
            filledArray = [filledArray, (inputArray(i-1)+1):(inputArray(i)-1)];
        end
            filledArray = [filledArray, inputArray(i)];
    end
end

function waveRho = updateRho(waveTime, wavePossible, wavePossibleRho)
    waveRho = [];
    for i = 1:length(waveTime)
        idx = find(wavePossible == waveTime(i), 1);
        if ~isempty(idx)
            waveRho = [waveRho, wavePossibleRho(idx)];
        else
            waveRho = [waveRho, NaN];
        end
    end
end

function [filteredArray, filteredRho] = removeShortSequences(inputArray, inputRho, minLength)
    filteredArray = []; filteredRho = [];
    currentSequence = []; currentRho = [];

    for i = 1:length(inputArray)
        if isempty(currentSequence) || inputArray(i) == currentSequence(end) + 1
            currentSequence = [currentSequence, inputArray(i)];
            currentRho = [currentRho, inputRho(i)];
        else
            if length(currentSequence) >= minLength
                filteredArray = [filteredArray, currentSequence];
                filteredRho = [filteredRho, currentRho];
            end
            currentSequence = inputArray(i);
            currentRho = inputRho(i);
        end
    end
    
    if length(currentSequence) >= minLength
        filteredArray = [filteredArray, currentSequence];
        filteredRho = [filteredRho, currentRho];
    end
end
