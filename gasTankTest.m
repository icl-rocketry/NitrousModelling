clear
clc

tank = GeometricNitrousTank(25+273.15,0,1,(0.25*pi*(15e-2)^2),0.8);
tank

mNitrousAdd = 1.5;
EOfNitrousAdd = mNitrousAdd * NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,30+273.15,FluidProperty.VAPOR_QUALITY,0);
tank.addAmountOfNitrous(mNitrousAdd, EOfNitrousAdd);
tank