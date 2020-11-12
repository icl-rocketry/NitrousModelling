clear
clc

%Example usage of GeometricNitrousTank.m

CELSIUS_TO_K = 273.15; %Use a constant to help with code readability
%Create nitrous tank variable with temp of 30C, Liquid mass of 6Kg and
%volume of 0.05m^3
internalTankHeight = 0.8; %Metres, very approximate geometry
internalTankCrossSectionA = 0.25 * pi * (150e-3)^2; %M^2, very approximate geometry
tank = GeometricNitrousTank(30+CELSIUS_TO_K,6,internalTankHeight,internalTankCrossSectionA,0.65);
%Maybe we want to show the tank pressure?
disp("Current tank pressure: "+tank.pressureAtBaseOfTank+" Pascal");
disp("Current tank mLiquid: "+tank.mLiquid);
disp("Start tank mTotal: "+tank.mTotalNitrous);
%Let's pretend the tank warned up by 5000J of heat for some reason
tank.addHeat(5000);
disp("Current tank temp: "+(tank.temp-CELSIUS_TO_K)+"C");
%Let's show the tank pressure again (It will now have changed)
disp("Current tank pressure: "+tank.pressureAtBaseOfTank+" Pascal");
disp("Current tank mLiquid: "+tank.mLiquid);

%Let's pretend we want to use up some oxidiser from our tank
mToDrain = 0.1;
Q = mToDrain*NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_INTERNAL_ENERGY,FluidProperty.TEMPERATURE,tank.temp,FluidProperty.PRESSURE,tank.pressureAtBaseOfTank);
tank.drainAmountOfLiquid(mToDrain); %0.1Kg
disp("Est required Q="+Q+"J");

%Values of tank pressure, temperature, liquid mass, etc... will now have changed.
%Let's show the user the remaining liquid mass and it's temperature
disp("Remaining liquid mass: "+tank.mLiquid+"Kg");
disp("Temperature of the tank: "+(tank.temp-CELSIUS_TO_K)+"C");
