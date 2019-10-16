%TO run you will need python 3.7 installed and to have executed the
%following commands (But replacing the path to python.exe)
%pyversion C:\Users\Edward\AppData\Local\Programs\Python\Python37\python.exe;
%[v,e] = pyversion; system([e,' -m pip install --user -U CoolProp']);

clear
clc
close all

%Change these parameters to alter the conditions of the chamber/injector
PDownstream = 30e5; %Pressure downstream of orifice (Chamber pressure)
AInjector1 = 0.25 * pi * (0.5e-3).^2; %Injector hole 1 cross section in m^3
AInjector2 = 0.25 * pi * (0.5e-3).^2; %Injector hole 2 cross section in m^3
%Discharge coefficient = "ratio of the actual discharge to the theoretical
%discharge". Ideally empirically determined
dischargeCoefficient1 = 0.8; %For injector orifice 1
dischargeCoefficient2 = 0.8; %For injector orifice 2

maxSimTime = 120; %120 sec

%Change these parameters to alter the conditions of the tank
initialInternalTankTemp = 27+273.15; %Starting tank temperature in kelvin
initialInternalNitrousMass = 6; %Kg
internalTankHeight = 0.8; %Metres, very approximate geometry
internalTankCrossSectionA = 0.25 * pi * (150e-3)^2; %M^2, very approximate geometry
characteristicPipeLength = 0.5; %0.5m characteristic pipe length

%Initialize our tank object and display
tank = GeometricNitrousTank(initialInternalTankTemp,initialInternalNitrousMass,internalTankHeight,internalTankCrossSectionA,0.65);
disp("Initial tank pressure: "+tank.vapourPressure);
disp("Initial tank liquid fraction: "+(tank.liquidHeight/tank.tankHeight));

disp("Simming drain in 0.25sec internals...");
i=0;
dt = 0.25; %sec
t(1) = 0;
while(tank.mLiquid > 0.01 && tank.vapourPressure > PDownstream && (i<1 || t(i) < maxSimTime))
    i=i+1;
    PTank(i) = tank.vapourPressure;
    t(i) = dt.*(i-1);
    if(tank.vapourPressure <= PDownstream)
        mdot1(i) = 0;
        mdot2(i) = 0;
        continue;
    end
    [~,~,~,~,G] = SaturatedNitrous.getDownstreamSaturatedNHNEFlowCond(0,tank.temp,tank.vapourPressure,PDownstream,0,characteristicPipeLength);
    mdot1(i) = G .* AInjector1 .* dischargeCoefficient1;
    mdot2(i) = G .* AInjector2 .* dischargeCoefficient2;
    mdotTotal = mdot1(i)+mdot2(i);
    h = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.PRESSURE,tank.vapourPressure,FluidProperty.VAPOR_QUALITY,0);
    tank.changeNitrousMassEnergy(-mdotTotal.*dt,-h.*mdotTotal.*dt);
%     disp("mdot1: "+mdot1+"Kg/s");
%     disp("mdot2: "+mdot2+"Kg/s");
    disp("t="+t(i));
    drawnow;
end

figure();
hold on;
plot(t,mdot1);
plot(t,mdot2);
title("Injector mass flow vs t");
ylabel("Kg/s");
xlabel("Time (s)");
legend('Injector 1','Injector 2');
hold off;
figure();
plot(t,PTank);
ylabel("Pressure (Pa)");
xlabel("Time (s)");
title("Tank pressure vs t");
