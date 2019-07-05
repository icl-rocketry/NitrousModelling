%Fluid properties. First arg is the input to cool props as defined in http://www.coolprop.org/coolprop/HighLevelAPI.html#table-of-string-inputs-to-propssi-function
classdef FluidProperty
   enumeration
       DENSITY('D'),
       SPECIFIC_ENTHALPY('Hmass'),
       PRESSURE('P'),
       VAPOR_QUALITY('Q'),
       SPECIFIC_ENTROPY('Smass'),
       TEMPERATURE('T'),
       SPECIFIC_INTERNAL_ENERGY('Umass'),
       SPEED_OF_SOUND('speed_of_sound'),
       THERMAL_CONDUCTIVITY('conductivity'),
       IDEAL_CP('Cp0mass'),
       CV('Cvmass'),
       CP('Cpmass'),
       PHASE_INDEX('Phase'),
       PHASE_IDENTIFICATION_PARAMETER('PIP'),
       CRITICAL_DENSITY('rhocrit'),
       SPECIFIC_GIBBS_ENERGY('Gmass'),
       ISENTROPIC_EXPANSION_COEFF('isentropic_expansion_coefficient'),
       ISOBARIC_EXPANSION_COEFF('isobaric_expansion_coefficient'),
       MOLAR_MASS('molarmass'),
       CRITICAL_PRESSURE('pcrit'),
       PRANDTL_NUMBER('Prandtl'),
       CRITICAL_TEMPERATURE('Tcrit'),
       VISCOSITY('viscosity'),
       COMPRESSIBILITY_FACTOR('Z');
   end
   properties
      coolPropInput;
   end
   methods
       function obj = FluidProperty(coolPropInput)
            obj.coolPropInput = coolPropInput;
       end
       function val = getCoolPropInput(obj)
          val = obj.coolPropInput; 
       end
   end
end