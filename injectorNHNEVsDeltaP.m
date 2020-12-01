clear
clc
close all

%Change these parameters to alter the conditions of the chamber/injector
PDownstream = 30e5; %Pressure downstream of orifice (Chamber pressure)
AInjector1 = 0.25 * pi * (1.5e-3).^2; %Injector hole 1 cross section in m^3
AInjector2 = 0.25 * pi * (1.5e-3).^2; %Injector hole 2 cross section in m^3

%Discharge coefficient = "ratio of the actual discharge to the theoretical
%discharge". Ideally empirically determined. THIS IS RATIO OF REAL TO NHNE
%MODEL NOT RATIO OF REAL TO INCOMPRESSIBLE
dischargeCoefficient1 = 0.8; %For injector orifice 1
dischargeCoefficient2 = 0.8; %For injector orifice 2

dP = linspace(5e5, 40e5, 15);
PUpstream = PDownstream+dP;

characteristicLength = 5e-3; %Characteristic length used for NHNE
XUpstream = 0; %Assume fully saturated liquid upstream - big bad assumption
vUpstream = 0; %Assume stagnation ahead of injector - big bad assumption

for i=1:length(PUpstream)
    %Assume upstream flow is ON saturation line and FULLY liquid (bold
    %simplifying assumption)
    TUpstream = SaturatedNitrous.getSaturationTemperature(PUpstream(i));
    h1(i) = NitrousFluidCoolProp.getProperty(FluidProperty.SPECIFIC_ENTHALPY, FluidProperty.PRESSURE, PUpstream(i), FluidProperty.VAPOR_QUALITY, XUpstream);
    [~,~,v2(i),h2(i),G] = SaturatedNitrous.getDownstreamSaturatedNHNEFlowCond(XUpstream,TUpstream,PUpstream(i),PDownstream,vUpstream,characteristicLength);
    dh(i) = h2(i) - h1(i);
    mdot1(i) = G .* AInjector1 .* dischargeCoefficient1;
    mdot2(i) = G .* AInjector2 .* dischargeCoefficient2;
    mdotTotal(i) = mdot1(i) + mdot2(i);
end

figure();
plot(dP, mdotTotal);
title("NHNE mdot for upstream saturated stagnated liquid to downstream P="+PDownstream+"Pa");
xlabel("dP (Pa)");
ylabel("mdot");