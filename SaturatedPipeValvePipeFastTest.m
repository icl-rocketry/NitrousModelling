clear
clc

spvp = SaturatedPipeValvePipeFast('preBakedData/saturatedGasPipeValveFlowRates.mat');
tic;
mdot = spvp.getMassFlow(50e5,1e5,1);
toc;
disp(mdot);