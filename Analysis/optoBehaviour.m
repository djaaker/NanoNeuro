function optoBehaviour(IntanBehaviour,parameters)

%% Comparing Reaction Time 
optoRT = [];
noOptoRT = [];

for i=1:size(IntanBehaviour.cueHitTrace,2)
    if IntanBehaviour.cueHitTrace(i).opto == 1
        optoRT = [optoRT;IntanBehaviour.cueHitTrace(i).reactionTime];
    else
        noOptoRT= [noOptoRT;IntanBehaviour.cueHitTrace(i).reactionTime];
    end
end

ranksum(optoRT,noOptoRT)

data = optoRT;
% manually copy noOptoRT into data

figure,customBoxplot(data);

