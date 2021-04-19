fSample = 1;
nSamples = 1024;
t = (0: nSamples - 1) / fSample;
fExp = [0.3, 0.32];
expWave = exp(1i * 2 * pi * fExp(1) * t) + exp(1i * 2 * pi * fExp(2) * t);
point = 20: 10: 50;
nPoints = length(point);
pNoise = 0.2;
nRps = 1e2;

psdMusic = cell(nPoints, nRps);
psdMusicMean = cell(nPoints, 1);
psdMusicStd = cell(nPoints, 1);
for iPoint = 1: nPoints
    for iRp = 1: nRps
        noise = sqrt(pNoise / 2) * (randn(1, point(iPoint)) + 1i * randn(1, point(iPoint)));
        noisyExp = expWave(1: point(iPoint)) + noise;
        [~, cor] = corrmtx(noisyExp, 14, 'mod');
        [psdMusic{iPoint, iRp}, f] = pmusic(cor, 2, [], fSample);
    end
    psdMusicMean{iPoint} = mean(cell2mat(psdMusic(iPoint, :)), 2);
    psdMusicStd{iPoint} = std(cell2mat(psdMusic(iPoint, :)), [], 2);
end


figure;
for iPoint = 1: nPoints
    subplot(nPoints, 2, 2 * (iPoint - 1) + 1);
    for iRp = 1: nRps
        irPlot = plot(f, psdMusic{iPoint, iRp}, 'k', 'LineWidth', 2);
        hold on;
    end
    meanPlot = plot(f, psdMusicMean{iPoint}, 'g', 'LineWidth', 2);
    grid on; grid minor;
    legend([irPlot, meanPlot], {'Individual', 'Mean'});
    title(['PSD estimate using MUSIC: N = ', num2str(point(iPoint))]);
    xlabel('Normalised frequency');
    ylabel('Pseudospectrum');
    xlim([0.25 0.40]);
end

for iPoint = 1: length(point)
    subplot(length(point), 2, 2 * iPoint);
    for iRp = 1: nRps
        irPlot = plot(f, psdMusic{iPoint, iRp}, 'k', 'LineWidth', 2);
        hold on;
    end
    stdPlot = plot(f, psdMusicStd{iPoint}, 'b', 'LineWidth', 2);
    grid on; grid minor;
    legend([irPlot, stdPlot], {'Individual', 'Standard deviation'});
    title(['Standard deviation of the MUSIC estimate: N = ', num2str(point(iPoint))]);
    xlabel('Normalised frequency');
    ylabel('Pseudospectrum');
    xlim([0.25 0.40]);
end