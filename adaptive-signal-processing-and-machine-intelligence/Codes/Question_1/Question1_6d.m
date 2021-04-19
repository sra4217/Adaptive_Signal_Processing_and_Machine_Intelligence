addpath('./')
[YhatO, YO] = regval(Bols);
[YhatP, YP] = regval(BPCR);

eO = sum(sum((YO-YhatO).^2))/numel(YO);
eP = sum(sum((YP-YhatP).^2))/numel(YP);