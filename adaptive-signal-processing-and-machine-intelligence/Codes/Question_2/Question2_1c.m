nSamples = 1e3;
nRps = 1e2;
coefAr = [0.1 0.8];
orderAr = length(coefAr);
variance = 0.25;
delay = 1;
step = [0.05; 0.01];
nSteps = length(step);
leak = 0;
nTransients = 5e2;

arModel = arima('AR', coefAr, 'Variance', variance, 'Constant', 0);
arSignal = simulate(arModel, nSamples, 'NumPaths', nRps);
arSignal = arSignal';

mse = zeros(nSteps, nRps);
for iStep = 1: nSteps
    for iRp = 1: nRps
        signal = arSignal(iRp, :);
        [group] = preprocessing(signal, orderAr, delay);
        [~, ~, error] = lms(group, signal, step(iStep), leak);
        mse(iStep, iRp) = mean(error(nTransients + 1: end) .^ 2);
    end
end
emse = mean(mse - variance, 2);
misadj = emse / variance;

cov = [25/27, 25/54; 25/54, 25/27];
misadjApprox = step / 2 * trace(cov);
fprintf('Misadjustment: %.4f    %.4f\n', misadj);
fprintf('Approximated misadjustment: %.4f   %.4f\n', misadjApprox);