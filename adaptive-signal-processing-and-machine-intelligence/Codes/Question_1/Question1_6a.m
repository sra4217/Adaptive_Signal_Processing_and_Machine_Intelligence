load('./data/PCAPCR.mat');

svdX = svd(X);
svdXN = svd(Xnoise);

figure
stem(svdX, 'x', 'linewidth', 2, 'markersize', 10);
xlabel('Singular Value Index')
ylabel('Singular Value Magnitude')
title('Singular Values for signal X')

figure
stem(svdXN, 'x', 'linewidth', 2, 'markersize', 10);
xlabel('Singular Value Index')
ylabel('Singular Value Magnitude')
title('Singular Values for Xnoise')

figure
stem((svdX-svdXN).^2, 'x', 'linewidth', 2, 'markersize', 10);
xlabel('Singular Value Index')
ylabel('Magnitude')
title('Squared Error')