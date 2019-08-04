%Class to represent any flow restriction that can be modelled via a simple
%(perhaps variable) Cv - flow coefficient. Eg. a valve or an orifice
classdef (Abstract) FlowCoeffFlowRestriction < FlowRestriction
    methods (Abstract)
        coeff = getFlowCoefficient(obj)
    end
    
    methods 
        function mdot = getMassFlowForPressureChange(obj,dP,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
            if (fluidType == FluidType.NITROUS_GENERAL || fluidType == FluidType.NITROUS_LIQUID || fluidType == FluidType.NITROUS_GAS) && XUpstream ~= 1 && XUpstream ~= 0
               rho = SaturatedNitrous.getDensity(XUpstream,TUpstream,PUpstream); 
            else
                if(fluidType == FluidType.NITROUS_GENERAL)
                    if(XUpstream == 1)
                        rho = FluidType.NITROUS_GAS.getDensity(TUpstream,PUpstream);
                    elseif(XUpstream == 0)
                        rho = FluidType.NITROUS_LIQUID.getDensity(TUpstream,PUpstream);
                    end
                else
                    rho = fluidType.getDensity(TUpstream,PUpstream);
                end
            end
            mdot = FlowCoefficient.getMassFlowFromCoeff(obj.getFlowCoefficient(),dP,rho);
        end
        
        function dP = getPressureChangeForMassFlow(obj,mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
            if (fluidType == FluidType.NITROUS_GENERAL || fluidType == FluidType.NITROUS_LIQUID || fluidType == FluidType.NITROUS_GAS) && XUpstream ~= 1 && XUpstream ~= 0
               rho = SaturatedNitrous.getDensity(XUpstream,TUpstream,PUpstream); 
            else
                if(fluidType == FluidType.NITROUS_GENERAL)
                    if(XUpstream == 1)
                        rho = FluidType.NITROUS_GAS.getDensity(TUpstream,PUpstream); 
                    elseif(XUpstream == 0)
                        rho = FluidType.NITROUS_LIQUID.getDensity(TUpstream,PUpstream); 
                    end
                else
                    rho = fluidType.getDensity(TUpstream,PUpstream); 
                end
            end
            dP = FlowCoefficient.getPressureChangeFromMassFlow(obj.getFlowCoefficient(),mdot,rho);
        end
    end
end