classdef SaturatedPipeValvePipeFast
    properties
        data; %The pre baked data
    end
    properties  (Constant)
        DATA_SPACING = 100e3;
        MAX_P = 72e5;
        MIN_P = 1e5;
    end
    
    methods
        function obj = SaturatedPipeValvePipeFast(datafileName)
            obj.data = load(datafileName,'-mat');
            obj.data = obj.data.data; %Extract map from struct
        end
        
        function mdot = getMassFlow(obj,upstreamPressure,downstreamPressure,valveAmt)
            upstreamPVals = obj.getNearbyDataPoints(upstreamPressure,obj.MAX_P,obj.MIN_P);
            mdots = zeros(1,length(upstreamPVals));
            for i=1:length(upstreamPVals)
                downstreamPVals = obj.getNearbyDataPoints(downstreamPressure,min(obj.MAX_P,upstreamPVals(i)),obj.MIN_P);
                innerMDots = zeros(1,length(downstreamPVals));
                for j=1:length(downstreamPVals)
                    innerMDots(j) = obj.getFlowRateOfDatapoint(upstreamPVals(i),downstreamPVals(j),valveAmt);
                end
                mdots(i) = ppval(spline(downstreamPVals,innerMDots),downstreamPressure);
            end
            mdot = ppval(spline(upstreamPVals,mdots),upstreamPressure);
        end
    end
    methods(Access=private)
        function mdot = getFlowRateOfDatapoint(obj,upstreamPressure,downstreamPressure,valveAmt)
            mdot = polyval(obj.data([num2str(upstreamPressure),'|',num2str(downstreamPressure)]),valveAmt);
        end
        
        function PVals = getNearbyDataPoints(obj,P,maxP,minP)
            nearest = round(P/obj.DATA_SPACING) * obj.DATA_SPACING;
            if(nearest > maxP)
               nearest = maxP; 
            end
            if(nearest < minP)
               nearest = minP; 
            end
            PVals(1) = nearest;
            P2 = nearest + obj.DATA_SPACING;
            P3 = nearest - obj.DATA_SPACING;
            if(P2 > maxP)
                P2 = P3;
                P3 = nearest - 2.*obj.DATA_SPACING;
            end
            if(P3 < minP)
                P3 = nearest + 2.*obj.DATA_SPACING;
            end
            PVals(2) = P2;
            PVals(3) = P3;
        end
    end
end