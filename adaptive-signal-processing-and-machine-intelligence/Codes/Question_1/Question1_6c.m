Bols = (Xnoise'*Xnoise)\Xnoise'*Y;
Yols = Xnoise*Bols;

rank = 3;
[U,S,V] = svd(Xnoise);
U = U(:,1:rank);
S = S(1:rank,1:rank);
V = V(:,1:rank);
BPCR = V/S*U'*Y;
Ypcr = Xnoise*BPCR;

eOLS = sum(sum((Y-Yols).^2))/numel(Y);
ePCR = sum(sum((Y-Ypcr).^2))/numel(Y);

YolsT = Xtest*Bols;
YpcrT = Xtest*BPCR;

eOLST = sum(sum((Ytest-YolsT).^2))/numel(Ytest);
ePCRT = sum(sum((Ytest-YpcrT).^2))/numel(Ytest);