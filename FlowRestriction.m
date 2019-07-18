%Class to represent any flow restriction that can be modelled via a simple
%(perhaps variable) Cv - flow coefficient. Eg. a valve or an orifice
classdef (Abstract) FlowRestriction < handle
    methods (Abstract)
        [T,P,X,vDownstream] = getDownstreamTemperaturePressureFromMassFlow(obj,mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
        [T,mdot,X,vDownstream] = getDownstreamTemperatureMassFlowFromPressureChange(obj,dP,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
        mdot = getMassFlowForPressureChange(obj,dP,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
        dP = getPressureChangeForMassFlow(obj,mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
    end
    
    methods 
        
    end
end