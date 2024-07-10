function [angVel, ang, waveCurl] = rotationalCharacteristics(waveTime, cav, curlz, sourceCurrent)

    % if numel(waveTime) > 1
    omega = nan(8,8,length(waveTime)); waveCurl = omega;
    for i = 1:length(waveTime)
        time = waveTime(i);
        omega(:,:,i) = cav(:,:,time); 
    end
    
    angVel = mean(omega(sourceCurrent(1),sourceCurrent(2),:));
    ang = angVel * length(waveTime) * 180/pi;
    waveCurl = mean(curlz(sourceCurrent(1),sourceCurrent(2),time));
end
