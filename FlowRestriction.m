%Class to represent any flow restriction that can be modelled via a simple
%(perhaps variable) Cv - flow coefficient. Eg. a valve or an orifice
classdef (Abstract) FlowRestriction < handle
    methods (Abstract)
        coeff = getFlowCoefficient(obj)
        [T,P,X] = getDownstreamTemperaturePressureFromMassFlow(obj,mdot,fluidType,TUpstream,PUpstream,XUpstream)
        [T,mdot,X] = getDownstreamTemperatureMassFlowFromPressureChange(obj,dP,fluidType,TUpstream,PUpstream,XUpstream)
        mdot = getMassFlowForPressureChange(obj,dP,fluidType,TUpstream,PUpstream,XUpstream)
        dP = getPressureChangeForMassFlow(obj,mdot,fluidType,TUpstream,PUpstream,XUpstream)
    end
    
    methods 
        
    end
end