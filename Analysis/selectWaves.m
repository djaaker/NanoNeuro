function [WavesOut] = selectWaves(Waves,stEvalP,spEvalP)

WavesOut = Waves;

for i=1:size(Waves,2)
    posDel = find(Waves(i).evaluationPoints < stEvalP | Waves(i).evaluationPoints > spEvalP);
    WavesOut(i).wavePresent = WavesOut(i).wavePresent(stEvalP:spEvalP);
    WavesOut(i).waveStart = WavesOut(i).waveStart(stEvalP:spEvalP);
    WavesOut(i).evaluationPoints(posDel) = [];
    WavesOut(i).PGD = WavesOut(i).PGD(stEvalP:spEvalP);
    WavesOut(i).vx = WavesOut(i).vx(stEvalP:spEvalP);
    WavesOut(i).vy = WavesOut(i).vy(stEvalP:spEvalP);
    WavesOut(i).rho(posDel) = [];
    WavesOut(i).waveTime(posDel,:) = [];
    WavesOut(i).source(posDel,:) = [];
    WavesOut(i).nWaves = WavesOut(i).nWaves - length(posDel);
    WavesOut(i).speed(posDel) = [];
    WavesOut(i).waveDir(posDel) = [];
    WavesOut(i).wavelength(posDel) = [];
    WavesOut(i).waveDuration(posDel) = [];
    WavesOut(i).waveAmp(posDel) = [];
%     WavesOut(i).waveFreq(posDel) = [];
end

