classdef FluidPhase
   enumeration
       LIQUID('liquid'),
       GAS('gas'),
       TWO_PHASE('twophase'),
       SUPERCRITICAL_LIQUID('supercritical_liquid'),
       SUPERCRITICAL_GAS('supercritical_gas'),
       SUPERCRITICAL('supercritical'),
       AUTO('not_imposed');
   end
   properties
      coolPropInput;
   end
   methods
       function obj = FluidPhase(coolPropInput)
            obj.coolPropInput = coolPropInput;
       end
       function val = getCoolPropInput(obj)
          val = obj.coolPropInput; 
       end
   end
end