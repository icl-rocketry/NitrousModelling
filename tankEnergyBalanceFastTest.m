% clear
% clc

tankBal = TankEnergyBalanceFast.INSTANCE;
% tankBal = TankEnergyBalanceFast('preBakedData/nitrousTankInternalEnergiesv3.mat');
mTest = 17.7 / 0.01767;
tic;
tankBal.getTemp(mTest+0.1,99500)
toc;
tic;
tankBal.getTemp(mTest-0.1,99500)
toc;