%Prefer to use simscape/simulink model
clear
clc
close all

pipe = FluidPipe(0.25*pi*(10e-3).^2,1);
PUpstream = 70e5;
TUpstreams = SaturatedNitrous.getSaturationTemperature(PUpstream);
XUpstream = 0;
vUpstream = 0;
PDownstream = 6995774:-0.5:6995770;
for j=1:length(TUpstreams)
    TUpstream = TUpstreams(j);
    for i=1:length(PDownstream)
        dP(i) = PDownstream(i) - PUpstream;
        [~,mdotOriginal(j,i),~,~] = pipe.getDownstreamTemperatureMassFlowFromPressureChange(dP(i),FluidType.NITROUS_GENERAL,TUpstream,PUpstream,XUpstream,vUpstream);
        disp("Approximate about "+dP(1));
        [T(i),mdot(j,i),X(i),vDownstream(i)] = pipe.getDownstreamTemperatureMassFlowFromDPAboutPt(dP(1),dP(i),FluidType.NITROUS_GENERAL,TUpstream,PUpstream,XUpstream,vUpstream);
        k(j,i) = mdot(j,i) / sqrt(abs(dP(i)));
        specificGravity = SaturatedNitrous.getDensity(X(i),T(i),PDownstream(i)) ./ FlowCoefficient.RHO_WATER;
        Q = ((mdot(j,i).*3600) / (SaturatedNitrous.getDensity(X(i),T(i),PDownstream(i)))); %m^3/hour
        dpBar = abs(dP(i)) / 1e5;
        kv(j,i) = Q * sqrt(specificGravity / dpBar);
    end
end

figure();
plot(sqrt(abs(dP)),mdot);
title('Sqrt pressure change vs mass flow (about pt)');
xlabel('Sqrt pressure change');
ylabel('Mass flow');

figure();
plot(sqrt(abs(dP)),mdotOriginal);
title('Sqrt pressure change vs mass flow (direct)');
xlabel('Sqrt pressure change');
ylabel('Mass flow');

figure();
plot(sqrt(abs(dP)),T);
title('Sqrt pressure change vs temperature');
xlabel('Sqrt pressure change');
ylabel('temperature');

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

for i=1:length(TUpstreams)
    [Tsats(i),PSats(i)] = SaturatedNitrous.getWhereIsentropicIntersectWithSaturationCurve(TUpstreams(i),PUpstream);
    dpSat(i,:) = PSats(i)-PDownstream(1,:);
end

figure();
hold on;
for j=1:length(TUpstreams)
    plot(PDownstream,kv(j,:));
end
hold off;
title('Downstream pressure vs metric flow coefficient');
xlabel('Downstream pressure');
ylabel('Metric flow coefficient');

figure();
hold on;
for j=1:length(TUpstreams)
    plot(mdot(j,:),kv(j,:));
end
hold off;
title('Mass flow vs metric flow coefficient');
xlabel('Mass flow');
ylabel('Metric flow coefficient');

figure();
hold on;
for j=1:length(TUpstreams)
    plot(log(mdot(j,:)),log(kv(j,:)));
end
hold off;
title('LOG Mass flow vs LOG metric flow coefficient');
xlabel('Log Mass flow');
ylabel('Log Metric flow coefficient');

figure();
hold on;
for j=1:length(TUpstreams)
    plot(dpSat(j,:),kv(j,:));
end
hold off;
title('SHIFTED Downstream pressure vs metric flow coefficient');
xlabel('SHIFTED Downstream pressure');
ylabel('Metric flow coefficient');

figure();
hold on;
for j=1:length(TUpstreams)
    plot(log(PDownstream),log(kv(j,:)));
end
hold off;
title('LOG Downstream pressure vs LOG metric flow coefficient');
xlabel('log Downstream pressure');
ylabel('log Metric flow coefficient');
% 
% figure();
% plot(PDownstream,k);
% title('Downstream pressure vs raw flow coefficient');
% xlabel('Downstream pressure');
% ylabel('raw flow coefficient');

% tic;
% dPForMassFlow = pipe.getPressureChangeForMassFlow(3,FluidType.NITROUS_GENERAL,TUpstream,PUpstream,XUpstream,vUpstream);
% toc;
% disp("For 3kg/s mass flow need dp of "+dPForMassFlow);
