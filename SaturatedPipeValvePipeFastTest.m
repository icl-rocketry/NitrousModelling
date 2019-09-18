clear
clc

spvp = SaturatedPipeValvePipeFast('preBakedData/saturatedGasPipeValveFlowRates.mat');
tic;
valvePos = 0:0.0025:1;
for i=1:length(valvePos)
    mdot(i) = spvp.getMassFlow(50e5,1e5,valvePos(i));
end
valvePosPolynomialFitCoeffs = polyfit(mdot,valvePos,20);
valvePosPolynomialFitCoeffs(length(valvePosPolynomialFitCoeffs)) = 0;
toc;
plot(valvePos,mdot);
inverseFunc = @(x) polyval(valvePosPolynomialFitCoeffs,x);
for j=1:length(mdot)
    valvePosCalc(j) = inverseFunc(mdot(j));
end
hold on;
plot(valvePosCalc,mdot);
legend('SPVPF','Inverted poly');
hold off;
tic;
disp(inverseFunc(0));
toc;
inverseFunc(spvp.getMassFlow(50e5,1e5,0.5))
%disp(mdot);