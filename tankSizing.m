%Script to aide finding tank geometry.
%Depends on NitrousTank.m and SaturatedNitrous.m
%By Eddie Brown

clear
clc
CELSIUS_TO_K = 273.15; %Use a constant to help with code readability

%Define tank geometry (Edit these variables to change your tank)
tankRadius = 0.065; %Radius of tank internal geometry in metres
tankHeight = 0.8;%1.35;%1.15; %Height of tank internal geometry in metres
massToDrain = 4; %The amount of liquid oxidiser required for the flight (Kg)
tankTemp = 30+CELSIUS_TO_K; %The tank starting temp in K
ventHeight = 0.75*tankHeight; %Max liquid height in tank, aka bottom of vent hole, in metres
disp("Vent height from top: "+(tankHeight-ventHeight)+"m");

%Calculate from the tank definition the variables required to define the
%tank for modelling it in flight
tankCrossSectionA = pi * tankRadius^2;
tankVol = tankCrossSectionA * tankHeight;
liquidHeight = ventHeight;
liquidVol = tankCrossSectionA * liquidHeight;
liquidDensity = SaturatedNitrous.getLiquidDensity(tankTemp);
liquidMass = liquidVol * liquidDensity;
disp("Tank volume: "+tankVol+"m^3"); %Display tank's volume
disp("Initial liquid Mass: "+liquidMass+"Kg"); %Display tank's initial liquid mass for this temperature

%Define the flight tank model
tank = NitrousTank(tankTemp,liquidMass,tankVol); 
disp("Initial vapour mass: "+tank.mVapour+"Kg"); %Display initial vapour mass
disp("Initial tank pressure: "+tank.pressure+" Pascal");
drawnow;

if(true)
    return;
end

%Drain from the tank the amount used throughout the flight
tank.drainAmountOfLiquid(massToDrain);

%Output the final state of the tank
disp("Remaining liquid mass: "+tank.mLiquid+"Kg"); %Display how much liquid mass is remaining
disp("End Temperature of the tank: "+(tank.temp-CELSIUS_TO_K)+"C"); %Display the end temp of the tank
disp("End pressure of the tank: "+tank.pressure+" Pascal"); %Display the end pressure of the tank
extraMassMarginPercent = (tank.mLiquid / massToDrain)*100; %Calculate how much extra mass is left compared to how much mass was taken from the tank
disp("Extra mass: "+extraMassMarginPercent+"%"); %Display the extra mass %

