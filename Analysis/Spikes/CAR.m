function [Spikes] = CAR(Spikes)
%common average referencing 

    Spikes.commAvgRef = mean(Spikes.hpSignal,1,"omitnan");
    Spikes.hpSpikes = Spikes.hpSignal - Spikes.commAvgRef;
end

