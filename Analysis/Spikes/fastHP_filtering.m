function Spikes = fastHP_filtering(intanData,Fs)
tic
global useGPU
if nargin < 2 || strcmp(Fs,'')
    Fs = 20000;
    disp(['Sampling rate set at ' num2str(Fs) ' Hz for filtering']);
end

Spike.Fs = Fs;
channel_num = size(intanData,1);
%High Pass Filtering
Fc = [300 3000];
Wn = Fc./(Fs/2);
b = fir1(5000,Wn,'bandpass');

%% GPU case
if useGPU
    buff = gpuArray(double(intanData)'); 
    datr2 = filter(b,1,buff); % causal forward filter
    datr2  = flipud(datr2); % reverse time
    datr2  = filter(b, 1, datr2); % causal forward filter again
    spikes = flipud(datr2); % reverse time back
    
    %Output GPU data
    Spikes.hpSpikes = gather(spikes);
else
    Spikes.hpSpikes = filtfilt(b,1,double(intanData)');
end
Spikes.hpSpikes = Spikes.hpSpikes';
Spikes.channel_num = channel_num;
toc
