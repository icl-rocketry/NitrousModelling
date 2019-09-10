clear
clc

CELSIUS_TO_K = 273.15; %Use a constant to help with code readability
tankRadius = 0.065; %Radius of tank internal geometry in metres
tankCrossSectionA = pi * tankRadius^2;
tankHeight = 1;
ambientTemp = CELSIUS_TO_K + 40;

tank = GeometricNitrousAndAirTank(30+CELSIUS_TO_K,7,0,tankHeight,tankCrossSectionA,0.8);
T1 = tank.temp;
disp("Tank temp: "+tank.temp);
disp("Tank pressure: "+tank.vapourPressure);
mFill = 50e-3;
disp("Liq drain rate: "+mFill+" Kg/s");
targetMass = 6;
uBefore = tank.getInternalEnergy()/tank.mTotalNitrous;
Q = tank.findHeatRateInputToKeepTempForDrainingLiquidRate(mFill);
tank.drainSmallAmountOfLiquid(mFill);
% Q = mOut * NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_INTERNAL_ENERGY,FluidProperty.TEMPERATURE,tank.temp,FluidProperty.PRESSURE,tank.getPressureAtHeight(0));
% disp("Estimated Q: "+Q);
dT = tank.temp - T1;
disp("dT 1 sec: "+dT);
% Q = tank.findHeatToMakeTemp(T1);
% disp("Estimated Q: "+Q);
disp("Tank temp: "+tank.temp);
tank.addHeat(Q);
disp("Tank temp: "+tank.temp);

dEdt = tank.findIntEnergyChangeRateForConstTemperatureWithFillRate(mFill);
disp("dEdt for constant T: "+dEdt+" W");
disp("Q rate req: "+Q+" W");
disp("----------------");

heatTransferCoeff = 0.2*205;
internalTankSurfaceArea = 2*tankCrossSectionA + pi*2*tankRadius*tankHeight;
QRateIntoTank = internalTankSurfaceArea .* heatTransferCoeff .* (ambientTemp - tank.temp);

mdotOut = calcMDotOutReq(QRateIntoTank,mFill,tank,tank.temp+2);
mdotIn = mdotOut+mFill;
disp("Req vent rate: "+mdotOut*1000+" g/sec");
disp("Req flow between tanks: "+mdotIn*1000+" g/sec");
disp("Fill efficiency (% of flow stays in tank): "+((100*(mFill))/(mdotIn))+"%");
timeToFill = targetMass / mFill;
disp("Time to fill: "+timeToFill+" sec ("+(timeToFill/60)+" mins)");
wastedNitrous = timeToFill*mdotOut;
disp("Wasted nitrous: "+wastedNitrous+" Kg");
disp("Used nitrous: "+(mdotIn*timeToFill)+" Kg");
disp("Overall nitrous efficiency: "+((100*targetMass)/(mdotIn*timeToFill))+" %");

function mdotOut = calcMDotOutReq(HeatTransferIn,mdotFill,tank,extTankT)
    hFromExt = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,extTankT,FluidProperty.VAPOR_QUALITY,0);
    hInt = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,tank.temp,FluidProperty.VAPOR_QUALITY,1);
    disp("hInt: "+hInt);
    disp("hFromExt: "+hFromExt);
    dEdt = tank.findIntEnergyChangeRateForConstTemperatureWithFillRate(mdotFill);
    disp("dE/dt: "+dEdt);
    mdotOut = (HeatTransferIn - dEdt + mdotFill.*hFromExt)/(hInt - hFromExt);
end
