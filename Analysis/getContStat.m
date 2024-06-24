function [p,h] = getContStat(x1,x2)

assert(size(x1,2)==size(x2,2), 'Size of the two inputs are not same')

for i=1:size(x1,2)
    [p(i), h(i)] = ranksum(x1(:,i), x2(:,i));
end

