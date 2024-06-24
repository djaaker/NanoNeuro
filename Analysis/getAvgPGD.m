function [avgPGD,avgPGDNull] = getAvgPGD(xgp,Waves,behaviourTrace,nullSpace,parameters)

if ~exist('nullSpace','var')
    nullSpace = 0;
end

nNullSpace = 10;

avgPGD= mean(vertcat(Waves.PGD),1);

if nullSpace == 1
    PGDShuffle = zeros(nNullSpace,size(xgp,2),size(xgp{1,1},3));
    for j=1:nNullSpace
        xgpshuffle = cellfun(@shuffle_channels, xgp, 'UniformOutput', false);
        PGDShuffle(j,:,:) = getPGD(xgpshuffle,behaviourTrace,parameters);
    end
    avgPGDNull = mean(PGDShuffle,[1 2]);
end
