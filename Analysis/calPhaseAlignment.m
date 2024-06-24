function [PA] = calPhaseAlignment(xgp,parameters)
% Ref: Spontaneous travelling cortical waves gate perception in behaving
% primates, Nature 2020

PA = zeros(size(xgp{1,1},1),size(xgp{1,1},2),size(xgp{1,1},3));

N = size(xgp,2);

xgpnorm = cellfun(@(s) s./abs(s), xgp, 'UniformOutput', false);

a = 0;

for t=1:size(xgp{1,1},3)
    for i=1:size(xgp{1,1},1)
        for j=1:size(xgp{1,1},2)
            for k=1:N
                a = a + xgpnorm{1,k}(i,j,t);
            end
            PA(i,j,t) = abs(sum(a)/N);
            a = 0;
        end
    end
end

