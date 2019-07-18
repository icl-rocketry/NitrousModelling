classdef FlowCoefficient
    properties(Constant)
        PASCAL_TO_PSI = 0.000145038;
        RHO_WATER = 997; %kg/m^3
        US_GALLONS_TO_CUBIC_METRES = 0.00378541;
    end
    methods (Static)  
        %All units are SI, but Cv is the imperial flow coefficient
        function m = getMassFlowFromCoeff(Cv,dP,rhoUpstream)
            dPPsi = -dP.*FlowCoefficient.PASCAL_TO_PSI;
            specificGravity = rhoUpstream ./ FlowCoefficient.RHO_WATER;
            Q = Cv .* sqrt(dPPsi./specificGravity); %In US gallons per min
            m = (Q.*FlowCoefficient.US_GALLONS_TO_CUBIC_METRES.*rhoUpstream)./60; %kg/sec
        end
        
        %All units are SI, but Cv is the imperial flow coefficient
        function dP = getPressureChangeFromMassFlow(Cv,mdot,rhoUpstream)
            specificGravity = rhoUpstream ./ FlowCoefficient.RHO_WATER;
            Q = ((mdot.*60) / (rhoUpstream)); %m^3/min
            Q = Q / FlowCoefficient.US_GALLONS_TO_CUBIC_METRES; %US gallons/min
            dPPsi = specificGravity .* (Q./Cv).^2;
            dP = -(dPPsi ./ FlowCoefficient.PASCAL_TO_PSI);
        end
    end
    
end