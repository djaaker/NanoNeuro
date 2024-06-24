function [PSDCh, f] = getAvgPSD(behaviourTrace,parameters)

PSDCh = zeros(size(behaviourTrace,2),parameters.rows*parameters.cols,(parameters.Fs/2)+1);

for trialno = 1:size(behaviourTrace,2)
    xf1 = behaviourTrace(trialno).rawLFP;
    for i=1:parameters.rows
        for j=1:parameters.cols
            a = squeeze(xf1(i,j,:));
            if sum(isnan(a))>0
                PSDCh(trialno,(i-1)*parameters.cols + j,:) = NaN;
            else
                [PSDCh(trialno,(i-1)*parameters.cols + j,:) ,f] = pwelch(a,600,0,1000,parameters.Fs);
            end
        end
    end
end

