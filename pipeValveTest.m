clear
clc
close all

upstreamP = 70e5;
upstreamT = SaturatedNitrous.getSaturationTemperature(upstreamP);
downstreamP = 1e5;
upstreamQuality = 0;
dischargeCoefficient = 0.8; %Ratio of actual mass flow to isentropic predicted
dP = downstreamP - upstreamP;

pipeInternalDiameter = 10.2108e-3;
pipe1 = FluidPipe(0.25*pi*(pipeInternalDiameter).^2,1);
valveOpenAmt = 0:0.025:1;
% valve = BallValve(12,valveOpenAmt(1));
% disp("Flow coeff: "+valve.getFlowCoefficient());
pipe2 = FluidPipe(0.25*pi*(pipeInternalDiameter).^2,1);
% pvp = PipeValvePipe(pipe1,valve,pipe2);

% valve = BallValve(0.2,0.025);
valve = BallValve(1.4,0.3);
pvp = PipeValvePipe(pipe1,valve,pipe2);
tic;
[TTest,mdotTest,XTest,vDownstreamTest] = pvp.getDownstreamTemperatureMassFlowFromPressureChange(dP,FluidType.NITROUS_GENERAL,upstreamT,upstreamP,upstreamQuality,0);
mdotTestCorrected = dischargeCoefficient .* mdotTest;
disp((mdotTestCorrected*1000)+" g/sec");
toc;
if(XTest > 0.999)
    a = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.SPEED_OF_SOUND,FluidProperty.TEMPERATURE,TTest,FluidProperty.PRESSURE,downstreamP);
elseif(XTest < 0.0001)
    a = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPEED_OF_SOUND,FluidProperty.TEMPERATURE,TTest,FluidProperty.PRESSURE,downstreamP);
else
    a = NitrousFluidCoolProp.getProperty(FluidProperty.SPEED_OF_SOUND,FluidProperty.TEMPERATURE,TTest,FluidProperty.VAPOR_QUALITY,XTest);
end
MTest = vDownstreamTest / a;
if true
    return;
end

tic;
parfor i=1:length(valveOpenAmt)
    try
        valve = BallValve(0.2,valveOpenAmt(i));
        pvp = PipeValvePipe(pipe1,valve,pipe2);
        [T(i),mdot(i),X(i),vDownstream(i)] = pvp.getDownstreamTemperatureMassFlowFromPressureChange(dP,FluidType.NITROUS_GENERAL,upstreamT,upstreamP,upstreamQuality,0);
        mdotCorrected(i) = dischargeCoefficient .* mdot(i);
        if(X(i) > 0.999)
        	a = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.SPEED_OF_SOUND,FluidProperty.TEMPERATURE,T(i),FluidProperty.PRESSURE,downstreamP);
        elseif(X(i) < 0.0001)
            a = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPEED_OF_SOUND,FluidProperty.TEMPERATURE,T(i),FluidProperty.PRESSURE,downstreamP);
        else
            a = NitrousFluidCoolProp.getProperty(FluidProperty.SPEED_OF_SOUND,FluidProperty.TEMPERATURE,T(i),FluidProperty.VAPOR_QUALITY,X(i));
        end
        M(i) = vDownstream(i) / a;
    catch excep
        disp("valveOpenAmt: "+valveOpenAmt(i));
        drawnow;
        rethrow(excep);
    end
end
toc;

figure();
plot(valveOpenAmt,mdotCorrected);
title('Valve open amt vs mass flow');
xlabel('Valve open amt (0->1)');
ylabel('Mass flow');

figure();
plot(valveOpenAmt,M);
title('Valve open amt vs exit mach');
xlabel('Valve open amt (0->1)');
ylabel('Exit mach');
% 
% figure();
% plot(valveOpenAmt,T);
% title('Valve open amt vs temperature');
% xlabel('Valve open amt (0->1)');
% ylabel('Temperature');

mdotPolynomialFitCoeffs = polyfit(valveOpenAmt,mdotCorrected,10);
% TPolynomialFitCoeffs = polyfit(valveOpenAmt,T,10);
% vDownstreamFitCoeffs = polyfit(valveOpenAmt,vDownstream,10);
% xFitCoeefs = polyfit(valveOpenAmt,X,10);