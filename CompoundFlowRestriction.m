%Represents flow restrictions connected in series
classdef CompoundFlowRestriction < FlowRestriction
    properties
        flowRestrictions;
    end
    methods
        %Varargin are classes they implement Flow Restriction that make up
        %this compound flow restriction. The order of varargin should be
        %from upstream to downstream
        function obj = CompoundFlowRestriction(varargin)
            obj.flowRestrictions = varargin;
        end
        
        function [T,P,X,vDownstream] = getDownstreamTemperaturePressureFromMassFlow(obj,mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
            T = TUpstream;
            P = PUpstream;
            X = XUpstream;
            vDownstream = vUpstream;
            
            for i=1:length(obj.flowRestrictions)
               [T,P,X,vDownstream] = obj.flowRestrictions{i}.getDownstreamTemperaturePressureFromMassFlow(mdot,fluidType,T,P,X,vDownstream);
            end
        end
        
        function [T,mdot,X,vDownstream] = getDownstreamTemperatureMassFlowFromPressureChange(obj,dP,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
            if(dP == 0)
               T = TUpstream;
               mdot = 0;
               X = XUpstream;
               vDownstream = vUpstream;
               return;
            end
            mdot = 0.1;
            dPMeasured = Inf;
            %Use flow coefficient converging method - Estimate flow
            %coefficient, use to compute mdot and do until converge within
            %100 pascal
            while(abs(dPMeasured - dP) > 100)
                [T,P,X,vDownstream] = obj.getDownstreamTemperaturePressureFromMassFlow(mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream);
                dPMeasured = P - PUpstream;
                %dPMeasured = obj.getPressureChangeForMassFlow(mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream);
                k = mdot / sqrt(dPMeasured);
                mdot = k * sqrt(dP);
            end
%             disp("STOP ITER");
%             drawnow;
%             pause(5);
%             
%             mdot = abs(real(fzero(@(X) dpErr(X),0.1)));
%             [T,~,X,vDownstream] = getDownstreamTemperaturePressureFromMassFlow(mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream);
%             
%             function err = dpErr(mdotGuess)
%                 mdotGuess = abs(real(mdotGuess));
%                 dPActual = obj.getPressureChangeForMassFlow(mdotGuess,fluidType,TUpstream,PUpstream,XUpstream,vUpstream);
%                 disp(mdotGuess+" "+dPActual+" ("+dP+")");
%                 drawnow;
%                 err = dPActual - dP;
%             end
        end
        
        function mdot = getMassFlowForPressureChange(obj,dP,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
            [~,mdot,~,~] = obj.getDownstreamTemperatureMassFlowFromPressureChange(dP,fluidType,TUpstream,PUpstream,XUpstream,vUpstream);
        end
        
        function dP = getPressureChangeForMassFlow(obj,mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
            [~,P,~,~] = obj.getDownstreamTemperaturePressureFromMassFlow(mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream);
            dP = P-PUpstream;
        end
    end
end