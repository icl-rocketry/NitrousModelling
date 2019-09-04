clear
clc

CELSIUS_TO_K = 273.15; %Use a constant to help with code readability
tankRadius = 0.065; %Radius of tank internal geometry in metres
tankCrossSectionA = pi * tankRadius^2;

tank = GeometricNitrousAndAirTank(20+CELSIUS_TO_K,7,0,1,tankCrossSectionA,0.8);
T1 = tank.temp;
disp("Tank temp: "+tank.temp);
mOut = 20e-3;
uBefore = tank.getInternalEnergy()/tank.mTotalNitrous;
tank.drainSmallAmountOfLiquid(mOut);
% Q = mOut * NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_INTERNAL_ENERGY,FluidProperty.TEMPERATURE,tank.temp,FluidProperty.PRESSURE,tank.getPressureAtHeight(0));
% disp("Estimated Q: "+Q);
dT = tank.temp - T1;
disp("dT: "+dT);
Q = tank.findHeatToMakeTemp(T1);
disp("Estimated Q: "+Q);
disp("Tank temp: "+tank.temp);
tank.addHeat(Q);
disp("Tank temp: "+tank.temp);

dEdt = tank.findIntEnergyChangeRateForConstTemperatureWithFillRate(mOut);
disp("dEdt: "+dEdt);
