clear
clc
close all

%Ball valve test
testValve = BallValve(20,0); %Cv = 12
angle = 1:0.01:90;
openAmt = angle ./ 90;
PUpstream = 47.1e5; %Pa
dP = 1e5; %1 bar
%mdot = 0.01; %10 gram a sec flow rate
for i=1:length(openAmt)
    testValve.setOpenAmtFraction(openAmt(i));
    Cv(i) = testValve.getFlowCoefficient(); 
    [T(i),mdot(i),X(i),~] = testValve.getDownstreamTemperatureMassFlowFromPressureChange(-dP,FluidType.NITROUS_LIQUID,290,PUpstream,0,0);
end

plot(angle,Cv);
title('Cv vs angle');
xlabel('Angle open');
ylabel('Cv');
figure();
plot(angle,mdot);
title('Mass flow vs angle for given dP');
xlabel('Angle open');
ylabel('Mass flow');
figure();
plot(angle,X);
title('downstream quality vs angle for given dP');
xlabel('Angle open');
ylabel('Quality');
figure();
plot(angle,T);
title('downstream temperature vs angle for given dP');
xlabel('Angle open');
ylabel('Quality');
figure();
plot(angle,PUpstream-dP);
title('Downstream pressure vs angle for given dP');
xlabel('Angle open');
ylabel('P');
hold on;
for i=1:length(T)
   PSat(i) = SaturatedNitrous.getVapourPressure(T(i));
end
plot(angle,PSat);
legend('Downstream pressure','Saturation pressure at temperature');