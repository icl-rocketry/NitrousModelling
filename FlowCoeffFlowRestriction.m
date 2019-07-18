%Class to represent any flow restriction that can be modelled via a simple
%(perhaps variable) Cv - flow coefficient. Eg. a valve or an orifice
classdef (Abstract) FlowCoeffFlowRestriction < FlowRestriction
    methods (Abstract)
        getFlowCoefficient(obj)
    end
    
    methods 
        function mdot = getMassFlowForPressureChange(obj,dP,fluidType,TUpstream,PUpstream,XUpstream)
            if (fluidType == FluidType.NITROUS_GENERAL || fluidType == FluidType.NITROUS_LIQUID || fluidType == FluidType.NITROUS_GAS) && XUpstream ~= 1 && XUpstream ~= 0
               rho = SaturatedNitrous.getDensity(XUpstream,TUpstream,PUpstream); 
            else
               rho = fluidType.getDensity(TUpstream,PUpstream); 
            end
            mdot = FlowCoefficient.getMassFlowFromCoeff(obj.getFlowCoefficient(),dP,rho);
        end
        
        function dP = getPressureChangeForMassFlow(obj,mdot,fluidType,TUpstream,PUpstream,XUpstream)
            if (fluidType == FluidType.NITROUS_GENERAL || fluidType == FluidType.NITROUS_LIQUID || fluidType == FluidType.NITROUS_GAS) && XUpstream ~= 1 && XUpstream ~= 0
               rho = SaturatedNitrous.getDensity(XUpstream,TUpstream,PUpstream); 
            else
               rho = fluidType.getDensity(TUpstream,PUpstream); 
            end
            dP = FlowCoefficient.getPressureChangeFromMassFlow(obj.getFlowCoefficient(),mdot,rho);
        end
    end
end