function trialWaves = setTrialWaves(waveTime,sourceCurrent,waveRho,angVel,waveCurl,ang)

trialWaves = struct();

% populate the struct
trialTime = 1:3001; wavePresent = zeros(1, length(trialTime));

waveTime = waveTime(~cellfun(@isempty, waveTime));
allTimes = cell2mat(cellfun(@(x) x(:)', waveTime(~cellfun(@isempty, waveTime)), 'UniformOutput', false));
wavePresent(allTimes) = 1;


trialWaves.waveTime = waveTime;
trialWaves.wavePresent = wavePresent;
trialWaves.source =  sourceCurrent;
trialWaves.rho = waveRho;
trialWaves.angVel = angVel;
trialWaves.curl =  waveCurl;
trialWaves.ang = ang;

trialWaves.start = cellfun(@(x) x(1), waveTime(~cellfun(@isempty, waveTime)));
trialWaves.end = cellfun(@(x) x(end), waveTime(~cellfun(@isempty, waveTime)));
trialWaves.duration = trialWaves.end - trialWaves.start;

trialWaves.quadrant = ceil(trialWaves.start / 16);
trialWaves.cue = floor(trialWaves.start / 1500);

% 
% xph = xf(:,:,waveTime);
% [ft,signIF] = instantaneous_frequency( xph, parameters.Fs );
% if ~isnan(signIF)
%     fprintf('Actual Frequency! %d\n',ft)
%     ft = ft * signIF; trialWaves.iFreq = ft;
% else 
%     trialWaves.iFreq(:,wave_cnt) = 0;
% end
% 

end

