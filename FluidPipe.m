classdef FluidPipe < FlowRestriction
    properties
       crossSectionA; 
       pipeLength;
    end
    
    methods
        function obj = FluidPipe(crossSectionA,length) %Constructor for object
            obj.crossSectionA = crossSectionA;
            obj.pipeLength = length;
        end
        
        function [T,P,X,v] = getDownstreamTemperaturePressureFromMassFlow(obj,mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
            dP = obj.getPressureChangeForMassFlow(mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream);
            P = PUpstream+dP;
            [T,~,X,v] = obj.getDownstreamTemperatureMassFlowFromPressureChange(dP,fluidType,TUpstream,PUpstream,XUpstream,vUpstream);
        end
        
        function [T,mdot,X,v] = getDownstreamTemperatureMassFlowFromPressureChange(obj,dP,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
            if(~(fluidType == FluidType.NITROUS_GENERAL || fluidType == FluidType.NITROUS_LIQUID || fluidType == FluidType.NITROUS_GAS))
               error('Unsupported fluid type. Only supports nitrous right now'); 
            end
            [X,T,v,~,G] = SaturatedNitrous.getDownstreamSaturatedNHNEFlowCond(XUpstream,TUpstream,PUpstream,PUpstream+dP,vUpstream,obj.pipeLength);
            mdot = G*obj.crossSectionA;
        end
        
        function mdot = getMassFlowForPressureChange(obj,dP,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
            if(~(fluidType == FluidType.NITROUS_GENERAL || fluidType == FluidType.NITROUS_LIQUID || fluidType == FluidType.NITROUS_GAS))
               error('Unsupported fluid type. Only supports nitrous right now'); 
            end
            [~,~,~,~,G] = SaturatedNitrous.getDownstreamSaturatedNHNEFlowCond(XUpstream,TUpstream,PUpstream,PUpstream+dP,vUpstream,obj.pipeLength);
            mdot = G*obj.crossSectionA;
        end
        
        function dP = getPressureChangeForMassFlow(obj,mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
            %Compute some equivalent to flow coefficient to get an
            %estimate, then iterate until close enouth to the final value
            dpPredicted = -200e3;
            mdotDataPt = -20;
            while(abs(mdotDataPt-mdot) > 5e-5)
                [~,mdotDataPt,~,~] = obj.getDownstreamTemperatureMassFlowFromPressureChange(dpPredicted,fluidType,TUpstream,PUpstream,XUpstream,vUpstream);
                k = mdotDataPt / sqrt(abs(dpPredicted));
                dpPredicted = -((mdot.^2)/(k.^2));
            end
            dP = dpPredicted;
%             toc;
%             disp("DP Predicted: "+dpPredicted);
%             disp("Err in prediction: "+errFunc(dpPredicted));
            
%             tic;
%             dP = -abs(real(fzero(@(X) errFunc(X),dpPredicted)));
%             toc;
%             disp("DP actual: "+dP);
%             disp("Err in actual: "+errFunc(dP));

%             function err = errFunc(dP)
%                 dP = -abs(real(dP));
%                 [~,mdotAchieved,~,~] = obj.getDownstreamTemperatureMassFlowFromPressureChange(dP,fluidType,TUpstream,PUpstream,XUpstream,vUpstream);
%                 err = mdot - mdotAchieved;
%             end
        end
    end
end