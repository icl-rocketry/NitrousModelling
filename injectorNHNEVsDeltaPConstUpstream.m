clear
clc
close all

%Change these parameters to alter the conditions of the chamber/injector
AInjector1 = 5.*0.25 * pi * (1.47e-3).^2; %Injector hole 1 cross section in m^3
AInjector2 = 5.*0.25 * pi * (1.47e-3).^2; %Injector hole 2 cross section in m^3

%Discharge coefficient = "ratio of the actual discharge to the theoretical
%discharge". Ideally empirically determined. THIS IS RATIO OF REAL TO NHNE
%MODEL NOT RATIO OF REAL TO INCOMPRESSIBLE
dischargeCoefficient1 = 0.765; %For injector orifice 1
dischargeCoefficient2 = 0.765; %For injector orifice 2
dischargeCoefficient1HEM = 0.765; %For injector orifice 1
dischargeCoefficient2HEM = 0.765; %For injector orifice 2

dP = linspace(1e5, 35e5, 15);
PUpstream = 40e5;
PDownstream = PUpstream-dP; %Pressure downstream of orifice (Chamber pressure)

characteristicLength = 5e-3; %Characteristic length used for NHNE
XUpstream = 0; %Assume fully saturated liquid upstream - big bad assumption
vUpstream = 0; %Assume stagnation ahead of injector - big bad assumption

for i=1:length(dP)
    %Assume upstream flow is ON saturation line and FULLY liquid (bold
    %simplifying assumption)
    TUpstream = SaturatedNitrous.getSaturationTemperature(PUpstream);
    h1(i) = NitrousFluidCoolProp.getProperty(FluidProperty.SPECIFIC_ENTHALPY, FluidProperty.PRESSURE, PUpstream, FluidProperty.VAPOR_QUALITY, XUpstream);
    [~,~,v2(i),h2(i),G] = SaturatedNitrous.getDownstreamSaturatedNHNEFlowCond(XUpstream,TUpstream,PUpstream,PDownstream(i),vUpstream,characteristicLength);
    [~,~,v2Hem(i),h2Hem(i),GHem] = SaturatedNitrous.getDownstreamIsentropicSaturatedHEMFlowCond(XUpstream,TUpstream,PUpstream,PDownstream(i),vUpstream);
    [~,~,~,~,GInc] = SaturatedNitrous.getDownstreamIncompressibleMassFlowTemp(XUpstream,PUpstream,TUpstream,PDownstream(i),vUpstream,true);
    dh(i) = h2(i) - h1(i);
    mdot1(i) = G .* AInjector1 .* dischargeCoefficient1;
    mdot2(i) = G .* AInjector2 .* dischargeCoefficient2;
    mdot1HEM(i) = GHem .* AInjector1 .* dischargeCoefficient1HEM;
    mdot2HEM(i) = GHem .* AInjector2 .* dischargeCoefficient2HEM;
    mdot1Inc(i) = GInc .* AInjector1 .* dischargeCoefficient1;
    mdot2Inc(i) = GInc .* AInjector2 .* dischargeCoefficient2;
    mdotTotal(i) = mdot1(i) + mdot2(i);
    mdotTotalHEM(i) = mdot1HEM(i) + mdot2HEM(i);
    mdotTotalInc(i) = mdot1Inc(i) +  mdot2Inc(i);
end

figure();
plot(dP, mdotTotal);
title('NHNE');
%title("NHNE mdot for upstream saturated stagnated liquid to downstream P="+PDownstream+"Pa");
xlabel("dP (Pa)");
ylabel("mdot");

figure();
plot(dP, mdotTotalHEM);
title('HEM');
%title("HEM mdot for upstream saturated stagnated liquid to downstream P="+PDownstream+"Pa");
xlabel("dP (Pa)");
ylabel("mdot");

figure();
plot(dP, mdotTotalInc);
title('Incompressible');
%title("HEM mdot for upstream saturated stagnated liquid to downstream P="+PDownstream+"Pa");
xlabel("dP (Pa)");
ylabel("mdot");