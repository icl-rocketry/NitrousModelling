clear
clc
close all

pipe = FluidPipe(0.25*pi*(10e-3).^2,1);
PUpstream = 60e5;
TUpstream = 295;
XUpstream = 0;
vUpstream = 0;
PDownstream = 60e5:-0.5e5:40e5;
for i=1:length(PDownstream)
    dP(i) = PDownstream(i) - PUpstream;
    [T(i),mdot(i),X(i),vDownstream(i)] = pipe.getDownstreamTemperatureMassFlowFromPressureChange(dP(i),FluidType.NITROUS_GENERAL,TUpstream,PUpstream,XUpstream,vUpstream); 
end

figure();
plot(sqrt(abs(dP)),mdot);
title('Sqrt pressure change vs mass flow');
xlabel('Sqrt pressure change');
ylabel('Mass flow');

figure();
plot(PDownstream,mdot);
title('Downstream pressure vs mass flow');
xlabel('Downstream pressure');
ylabel('Mass flow');

figure();
plot(PDownstream,T);
title('Downstream pressure vs Temperature');
xlabel('Downstream pressure');
ylabel('Temperature');

figure();
plot(PDownstream,X);
title('Downstream pressure vs downstream quality');
xlabel('Downstream pressure');
ylabel('Vapour Quality');

figure();
plot(PDownstream,vDownstream);
title('Downstream pressure vs vel');
xlabel('Downstream pressure');
ylabel('Downstream vel');

tic;
dPForMassFlow = pipe.getPressureChangeForMassFlow(3,FluidType.NITROUS_GENERAL,TUpstream,PUpstream,XUpstream,vUpstream);
toc;
disp("For 3kg/s mass flow need dp of "+dPForMassFlow);
