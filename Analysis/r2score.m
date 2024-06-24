function [r2score] = r2score(y1,y2)

assert((size(y1,2) == 1) && (size(y2,2) == 1),'Pass 1-D data structures')
assert(size(y1,1) == size(y2,1), 'Both data sets should be same length')

y1mean = mean(y1,1);

SSres = sum((y1-y2).^2,'all');
SStot = sum((y1-y1mean).^2,'all');

r2score = 1-(SSres/SStot);

end

