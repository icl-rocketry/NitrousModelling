clear
clc
close all

upstreamP = 24e5;
upstreamT = SaturatedNitrous.getSaturationTemperature(upstreamP);
downstreamP = 7e5;
upstreamQuality = 1;
dP = downstreamP - upstreamP;

pipe1 = FluidPipe(0.25*pi*(10e-3).^2,1);
valveOpenAmt = 0:0.025:1;
% valve = BallValve(12,valveOpenAmt(1));
% disp("Flow coeff: "+valve.getFlowCoefficient());
pipe2 = FluidPipe(0.25*pi*(10e-3).^2,1);
% pvp = PipeValvePipe(pipe1,valve,pipe2);

tic;
parfor i=1:length(valveOpenAmt)
    try
        valve = BallValve(12,valveOpenAmt(i));
        pvp = PipeValvePipe(pipe1,valve,pipe2);
        [T(i),mdot(i),X(i),vDownstream(i)] = pvp.getDownstreamTemperatureMassFlowFromPressureChange(dP,FluidType.NITROUS_GENERAL,upstreamT,upstreamP,upstreamQuality,0);
    catch excep
        disp("valveOpenAmt: "+valveOpenAmt(i));
        drawnow;
        rethrow(excep);
    end
end
toc;

figure();
plot(valveOpenAmt,mdot);
title('Valve open amt vs mass flow');
xlabel('Valve open amt (0->1)');
ylabel('Mass flow');
% 
% figure();
% plot(valveOpenAmt,T);
% title('Valve open amt vs temperature');
% xlabel('Valve open amt (0->1)');
% ylabel('Temperature');

mdotPolynomialFitCoeffs = polyfit(valveOpenAmt,mdot,10);
% TPolynomialFitCoeffs = polyfit(valveOpenAmt,T,10);
% vDownstreamFitCoeffs = polyfit(valveOpenAmt,vDownstream,10);
% xFitCoeefs = polyfit(valveOpenAmt,X,10);