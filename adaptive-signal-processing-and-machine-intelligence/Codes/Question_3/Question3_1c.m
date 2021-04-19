fSample = 1e4;
fPower = 50;
nSamples = 1e3;
nPhases = 3;
amplitude = ones(nPhases, 1);
phaseShift = [0; -2/3 * pi; 2/3 * pi];
phaseInit = 0;
t = (0: nSamples - 1) / fSample;

balancedABC = amplitude .* cos(2 * pi * fPower * t + phaseInit + phaseShift);
balancedZeroAlphaBeta = clarke(balancedABC);
balancedClarke = balancedZeroAlphaBeta(2, :) + 1i * balancedZeroAlphaBeta(3, :);
[circularityBalanced, ~] = circularity(balancedClarke);

ampDif = 0.2: 0.2: 0.8;
nAmps = length(ampDif);
unbalancedClarkeAmp = cell(nAmps, 1);
circularityUnbalancedAmp = zeros(nAmps, 1);
for iAmp = 1: nAmps
    unbalancedABCAmp = (amplitude + [-ampDif(iAmp); 0; ampDif(iAmp)]) .* cos(2 * pi * fPower * t + phaseInit + phaseShift);
    unbalancedZeroAlphaBetaAmp = clarke(unbalancedABCAmp);
    unbalancedClarkeAmp{iAmp} = unbalancedZeroAlphaBetaAmp(2, :) + 1i * unbalancedZeroAlphaBetaAmp(3, :);
    [circularityUnbalancedAmp(iAmp), ~] = circularity(unbalancedClarkeAmp{iAmp});
end

phaseDif = 0.05: 0.05: 0.2;
nPhases = length(phaseDif);
unbalancedClarkePhase = cell(nPhases, 1);
circularityUnbalancedPhase = zeros(nPhases, 1);
for iPhase = 1: nPhases
    phaseDelay = [0; -phaseDif(iPhase) * pi; phaseDif(iPhase) * pi];
    unbalancedABCPhase = amplitude .* cos(2 * pi * fPower * t + phaseInit + phaseShift + phaseDelay);
    unbalancedZeroAlphaBetaPhase = clarke(unbalancedABCPhase);
    unbalancedClarkePhase{iPhase} = unbalancedZeroAlphaBetaPhase(2, :) + 1i * unbalancedZeroAlphaBetaPhase(3, :);
    [circularityUnbalancedPhase(iPhase), ~] = circularity(unbalancedClarkePhase{iPhase});
end

legendStr = cell(nAmps + 1, 1);
figure;
subplot(1, 2, 1);
scatter(real(balancedClarke), imag(balancedClarke), 'k');
legendStr{1} = sprintf('Balanced \\rho = %.2f', circularityBalanced);
hold on;
for iAmp = 1: nAmps
    scatter(real(unbalancedClarkeAmp{iAmp}), imag(unbalancedClarkeAmp{iAmp}));
    legendStr{iAmp + 1} = sprintf('\\DeltaV = %.1f \\rho = %.2f', ampDif(iAmp), circularityUnbalancedAmp(iAmp));
    hold on;
end
legend(legendStr);
title('Circularity diagram with unbalanced magnitudes');
xlabel('Real part');
ylabel('Imaginary part');
xlim([-2 2]);
ylim([-2 2]);
set(gcf, 'position', [10, 10, 500, 500])
legendStr = cell(nPhases + 1, 1);
subplot(1, 2, 2);
scatter(real(balancedClarke), imag(balancedClarke), 'k');
legendStr{1} = sprintf('Balanced \\rho = %.2f', circularityBalanced);
hold on;
for iPhase = 1: nPhases
    scatter(real(unbalancedClarkePhase{iPhase}), imag(unbalancedClarkePhase{iPhase}));
    legendStr{iPhase + 1} = sprintf('\\Delta\\phi = %.2f\\pi \\rho = %.2f', phaseDif(iPhase), circularityUnbalancedPhase(iPhase));
    hold on;
end
legend(legendStr);
title('Circularity diagram with unbalanced phases');
xlabel('Real part');
ylabel('Imaginary part');
xlim([-2 2]);
ylim([-2 2]);
set(gcf, 'position', [10, 10, 500, 500])