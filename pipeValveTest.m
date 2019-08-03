clear
clc
close all

pipe1 = FluidPipe(0.25*pi*(10e-3).^2,1);
valveOpenAmt = 0:0.0005:1;
% valve = BallValve(12,valveOpenAmt(1));
% disp("Flow coeff: "+valve.getFlowCoefficient());
pipe2 = FluidPipe(0.25*pi*(10e-3).^2,1);
% pvp = PipeValvePipe(pipe1,valve,pipe2);
tic;
parfor i=1:length(valveOpenAmt)
    valve = BallValve(12,valveOpenAmt(i));
    pvp = PipeValvePipe(pipe1,valve,pipe2);
    [T(i),mdot(i),X(i),vDownstream(i)] = pvp.getDownstreamTemperatureMassFlowFromPressureChange(-40e5,FluidType.NITROUS_GENERAL,303.15,6.32e6,0,0);
end
toc;

figure();
plot(valveOpenAmt,mdot);
title('Valve open amt vs mass flow');
xlabel('Valve open amt (0->1)');
ylabel('Mass flow');