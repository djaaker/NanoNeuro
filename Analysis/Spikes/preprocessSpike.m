function output = preprocessSpike(data,Fs)
Fc = [500 3000];
Wn = Fc./(Fs/2);
b = fir1(5000,Wn,'bandpass');
rawspikeTrace = filtfilt(b,1,double(data)');
rawspikeTrace = rawspikeTrace';
commonModeAvg = rawspikeTrace-mean(rawspikeTrace,1,"omitnan");
whitenedSpikeTrace = commonModeAvg;
output.rawspikeTrace = rawspikeTrace;
output.whitenedSpikeTrace = whitenedSpikeTrace - mean(whitenedSpikeTrace,2);
end