%function pca_wave(behaviour,trial)

total_waves = 0; 
X = zeros(1,86);

for trial = 1:143
    trial = 1; waves = WavesBaseline.wavesHit(trial); 
    max_length = max(waves(trial).waveDuration); nWaves = waves(trial).nWaves;
    max_length = 40; 
    total_waves = total_waves + nWaves;
    
    pgd = zeros(nWaves, max_length);
    rho = zeros(nWaves, max_length);
    
    for j = 1:nWaves
        start = waves.waveTime(j,1); stop = waves.waveTime(j,2);
        pgd_dur = stop - start + 1; rho_dur = length(waves.rho{1,j});
    
        pgd(j,1:pgd_dur) = waves.PGD(start:stop);
        rho(j,1:rho_dur) = waves.rho{1,j};
    end
    
    wave_vars = nan(nWaves, 6);
    
    for j = 1:nWaves
        wave_vars(j,1) = waves.speed(j);
        %wave_vars(j,2) = waves.speed(j);
        wave_vars(j,2) = waves.source(j,1);
        wave_vars(j,3) = waves.waveDir(j);
        wave_vars(j,4) = waves.wavelength(j);
        wave_vars(j,5) = waves.waveDuration(j);
        wave_vars(j,6) = waves.waveAmp(j);
    end
    
    X_trial = [rho pgd wave_vars];
    X = [X; X_trial];
end

[coeff,score,latent] = pca(X);
cm = jet(500);

figure(1);
hold on 

for j = 1:500
    plot3(score(j,1),score(j,2),score(j,3),'.','Color',cm(j,:))
    j
end

xlabel('PC1')
ylabel('PC2')
zlabel('PC3')