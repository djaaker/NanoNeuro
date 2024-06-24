function [PGDfreqHit,PGDfreqMiss] = getPGDFreqBand(LFPdatacube,IntanBehaviour,nullSpace,parameters)

if ~exist('nullSpace','var')
    nullSpace = 0;
end

nNullSpace = 100;
nMax = 19;
PGDfreqHit = zeros(nMax+1,1);
PGDfreqMiss = zeros(nMax+1,1);

for i=1:nMax
    xf = bandpass_filter(LFPdatacube,i*5,(i+1)*5,4,1000);
    [xgp, ~] = generalized_phase(xf,1000,0);
    if nullSpace == 0
        PGDfreqHit1 = getPGD(xgp,IntanBehaviour.cueHitTrace,parameters);
        PGDfreqMiss1 = getPGD(xgp,IntanBehaviour.cueMissTrace,parameters);
        PGDfreqHit(i) = mean(PGDfreqHit1,'all');
        PGDfreqMiss(i) = mean(PGDfreqMiss1,'all');
    else
        PGDfreqHitShuffle = zeros(nNullSpace,IntanBehaviour.nCueHit,size(IntanBehaviour.cueHitTrace(2).trace,1));
        PGDfreqMissShuffle = zeros(nNullSpace,IntanBehaviour.nCueMiss,size(IntanBehaviour.cueMissTrace(2).trace,1));
        for j=1:nNullSpace
            PGDfreqHitShuffle(j,:,:) = getPGD(shuffle_channels(xgp),IntanBehaviour.cueHitTrace,parameters);
            PGDfreqMissShuffle(j,:,:) = getPGD(shuffle_channels(xgp),IntanBehaviour.cueMissTrace,parameters);
        end
        PGDfreqHit(i) = mean(PGDfreqHitShuffle,'all');
        PGDfreqMiss(i) = mean(PGDfreqMissShuffle,'all');
    end
end

