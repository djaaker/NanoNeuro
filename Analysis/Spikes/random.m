%% Spikes
Spikes = fastHP_filtering(Intan.allIntan,20000);
Spikes = CAR(Spikes);
Spikes = findSpikes(Spikes,3.5);

spikeCh = 18;

figure();plot(Spikes.hpSpikes(spikeCh,:));
times = 1000*linspace(1,111,111)/Fs;
figure();plot(times,filtSpikeChData(8889340:8889450),'Color','r','LineWidth',2)
ylabel('Voltage (\muV)');
xlabel('Time (ms)');
box off

