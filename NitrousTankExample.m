clear
clc

%Example usage of NitrousTank.m, make sure SaturatedNitrous.m is in the
%same folder aswell or will not work

CELSIUS_TO_K = 273.15; %Use a constant to help with code readability
%Create nitrous tank variable with temp of 30C, Liquid mass of 6Kg and
%volume of 0.05m^3
tank = NitrousTank(30+CELSIUS_TO_K,6,0.05);
%Maybe we want to show the tank pressure?
disp("Current tank pressure: "+tank.pressure+" Pascal");
disp("Current tank mLiquid: "+tank.mLiquid);
disp("Start tank mTotal: "+tank.mTotal);
%Let's pretend the tank warned up by 1deg for some reason
tank.temp = tank.temp+1; %Increase temp of tank by 1
%Let's show the tank pressure again (It will now have changed)
disp("Current tank pressure: "+tank.pressure+" Pascal");
disp("Current tank mLiquid: "+tank.mLiquid);

%Let's pretend we want to use up some oxidiser from our tank - Note that
%the amount of liquid in the tank will decrease by more than this due to
%the interactions going on in the tank (Is modelled by the tank variable so
%don't worry). Behaviour if you try and drain more liquid from the tank
%than it contains isn't explicitly defined, but you'll probably get an
%error about the 'mass drainage solver' failing. Best to just check that
%you are not draining more than exists. Also worth noting that this works
%by iteratively solving for small mass changes dm that sum up to this total
%amount, so entering a large amount (Such as 4Kg) means it may take a few seconds to
%compute whereas a tiny amount will be very fast.
Q = 0.1*NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_INTERNAL_ENERGY,FluidProperty.TEMPERATURE,tank.temp,FluidProperty.PRESSURE,tank.pressure);
tank.drainAmountOfLiquid(0.1); %0.1Kg
disp("Est required Q="+Q);

%Values of tank pressure, temperature, liquid mass, etc... will now have changed.
%Let's show the user the remaining liquid mass and it's temperature
disp("Remaining liquid mass: "+tank.mLiquid+"Kg");
disp("Temperature of the tank: "+(tank.temp-CELSIUS_TO_K)+"C");
