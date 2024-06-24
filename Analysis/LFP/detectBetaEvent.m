function [BetaEvent] = detectBetaEvent(LFP,behaviourTrace,parameters)

% Reference - https://academic.oup.com/brain/article/140/11/2968/4430808#112576655

assert( ~isreal(LFP), 'pass hilbert transform of the LFP');

% Checking the structure of LFP pass. 
dimension = size(LFP);
if numel(dimension)==2
    meanLFP = mean(abs(LFP),1);
    disp('Shank data passed');
elseif numel(dimension)==3
    meanLFP = mean(abs(LFP),[1,2],'omitnan');
    disp('Grid data passed');
end

betaBurstThres = 2.5; % 2.5 * std dev 
betaAttenThres = 2.5; % -2.5 * std dev
minInterRippleInterval = 10 * parameters.Fs/1000; % in points
minBetaDuration = 30 * parameters.Fs/1000; % in points
maxBetaDuration = 210 * parameters.Fs/1000; % in points

%% Getting LFP for each trial and channels 
rows = size(LFP,1);
cols = size(LFP,2);
nlength = size(behaviourTrace(2).trace,2);
nTrials = size(behaviourTrace,2);

for i=1:nTrials
    BetaEvent(i).xamp = meanLFP(:,behaviourTrace(i).LFPIndex(1):behaviourTrace(i).LFPIndex(end));
    BetaEvent(i).xpow = BetaEvent(i).xamp.^2;
end

%% Normalizing and getting thresholds
normByTrial = 1; % 1 - normalize power by noise in each trial, 0 - normalize power by noise across all trials

if normByTrial == 1
    for i=1:nTrials
        BetaEvent(i).xnorm = (BetaEvent(i).xpow - mean(BetaEvent(i).xpow,'all','omitnan'))/std(BetaEvent(i).xpow,'omitnan');
    end
else
    meanPower = mean(horzcat(BetaEvent(1:end).xpow),'all','omitnan');
    stdPower = std(horzcat(BetaEvent(1:end).xpow),'omitnan');
    for i=1:nTrials
        BetaEvent(i).xnorm = (BetaEvent(i).xpow - meanPower)/stdPower;
    end
end

%% Getting beta burst segments 
for i=1:nTrials
    burstSegments = findThresSeg(BetaEvent(i).xnorm,betaBurstThres);
    attenSegments = findThresSeg(-1*BetaEvent(i).xnorm,betaAttenThres);
    BetaEvent(i).betaBurstPresent = zeros(1,size(BetaEvent(i).xnorm,2));
    % rejecting segments that are less then minimum time and more than
    % maximum time
    segmentSize = diff(burstSegments,1,2);
    rejectIndex1 = find(segmentSize<minBetaDuration); % minimum duration
    rejectIndex2 = find(segmentSize>maxBetaDuration); % maximum duration
    rejectIndex = [rejectIndex1;rejectIndex2];
    burstSegments(rejectIndex,:) = [];
    if isempty(burstSegments)
        disp(['No beta burst detect for Trial ' num2str(i)])
    else
        % Joining segments if minimum interval between segments is less
        % than minInterRippleInterval
        % This won't work if there are three continous segments with short
        % intervals. Need to modify 
        if size(burstSegments,1)>1
            intervals = diff(reshape(burstSegments',[],1));
            intervals = intervals(2:2:end);
            indx1 = find(intervals<minInterRippleInterval);
            for j=1:numel(indx1)
                burstSegments(indx1(j),2) = burstSegments(indx1(j)+1,2);
                burstSegments(indx1(j)+1,:) = [];
            end
            burstSegments = removeNaNRows(burstSegments);
        end
    end
    BetaEvent(i).bursts = burstSegments;
    for jj=1:size(burstSegments,1)
        BetaEvent(i).betaBurstPresent(burstSegments(jj,1):burstSegments(jj,2)) = 1;
    end

    segmentSize = diff(attenSegments,1,2);
    rejectIndex1 = find(segmentSize<minBetaDuration); % minimum duration
    rejectIndex2 = find(segmentSize>maxBetaDuration); % maximum duration
    rejectIndex = [rejectIndex1;rejectIndex2];
    attenSegments(rejectIndex,:) = [];
    if isempty(attenSegments)
        % disp(['No beta attenuations detect for Trial ' num2str(i)])
    else
        % Joining segments if minimum interval between segments is less
        % than minInterRippleInterval
        % This won't work if there are three continous segments with short
        % intervals. Need to modify 
        if size(attenSegments,1)>1
            intervals = diff(reshape(attenSegments',[],1));
            intervals = intervals(2:2:end);
            indx1 = find(intervals<minInterRippleInterval);
            for j=1:numel(indx1)
                attenSegments(indx1(j),2) = attenSegments(indx1(j)+1,2);
                attenSegments(indx1(j)+1,:) = [];
            end
            attenSegments = removeNaNRows(attenSegments);
        end
    end
    BetaEvent(i).attens = attenSegments;
    if parameters.opto == 1
        BetaEvent(i).opto = behaviourTrace(i).opto;
    end
end