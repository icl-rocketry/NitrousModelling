clear
clc
close all

%Ball valve test
testValve = BallValve(12,0); %12
angle = 1:0.01:90;
openAmt = angle ./ 90;
mdot = 0.01; %10 gram a sec flow rate
for i=1:length(openAmt)
    testValve.setOpenAmtFraction(openAmt(i));
    Cv(i) = testValve.getFlowCoefficient(); 
    dP(i) = testValve.getPressureChangeForMassFlow(mdot,FluidType.NITROUS_LIQUID,290,47.1e5,0,0);
    [T(i),P(i),X(i),~] = testValve.getDownstreamTemperaturePressureFromMassFlow(mdot,FluidType.NITROUS_LIQUID,290,47.1e5,0,0);
    %mdot(i) = testValve.getMassFlowForPressureChange(dP,FluidType.NITROUS_LIQUID,290,60e5,0);
end

plot(angle,Cv);
title('Cv vs angle');
xlabel('Angle open');
ylabel('Cv');
figure();
plot(angle,dP);
title('pressure drop vs angle for 10gram/sec flow');
xlabel('Angle open');
ylabel('pressure drop');
figure();
plot(angle,X);
title('downstream quality vs angle for 10gram/sec flow');
xlabel('Angle open');
ylabel('Quality');
figure();
plot(angle,T);
title('downstream temperature vs angle for 10gram/sec flow');
xlabel('Angle open');
ylabel('Quality');
figure();
plot(angle,P);
title('Downstream pressure vs angle for 10gram/sec flow');
xlabel('Angle open');
ylabel('P');
hold on;
for i=1:length(T)
   PSat(i) = SaturatedNitrous.getVapourPressure(T(i));
end
plot(angle,PSat);
legend('Downstream pressure','Saturation pressure at temperature');