%Function to generate a FULL nitrous tank for a given take off temp (K), tank Ullage (Fraction 0->1), tank height (m) and tank cross section area (m^2) 
function tank = tankGen(tempTakeoff,tankUllage,tankHeight,tankCrossSectionA)
    PVap = SaturatedNitrous.getVapourPressure(tempTakeoff); %Vapour pressure at this temp
    liqMass = FluidType.NITROUS_LIQUID.getDensity(tempTakeoff,PVap) * tankHeight * (1-tankUllage) * tankCrossSectionA;
    vapMass = FluidType.NITROUS_GAS.getDensity(tempTakeoff,PVap) * tankHeight * (tankUllage) * tankCrossSectionA;
    %Generate the tank model
    tank = GeometricNitrousTank(tempTakeoff,(liqMass+vapMass),tankHeight,tankCrossSectionA,(1-tankUllage)*tankHeight);
end