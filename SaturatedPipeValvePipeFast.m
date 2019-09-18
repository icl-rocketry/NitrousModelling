classdef SaturatedPipeValvePipeFast < matlab.mixin.Copyable
    properties
        data; %The pre baked data
        dataInverted;
        hasInvertedData = false;
    end
    properties  (Constant)
        DATA_SPACING = 100e3;
        MAX_P = 72e5;
        MIN_P = 1e5;
    end
    
    methods
        function obj = SaturatedPipeValvePipeFast(datafileName,inverseDataFileName)
            obj.data = load(datafileName,'-mat');
            obj.data = obj.data.data; %Extract map from struct
            if ~exist('inverseDataFileName','var')
               obj.hasInvertedData = false;
            else
               obj.dataInverted = load(inverseDataFileName,'-mat');
               obj.dataInverted = obj.dataInverted.data;
               obj.hasInvertedData = true;
            end
        end
        
        function valveAmt = getValvePos(obj,upstreamPressure,downstreamPressure,mdot)
%             disp("Finding valve pos for mdot of "+mdot);
            upstreamPVals = obj.getNearbyDataPoints(upstreamPressure,obj.MAX_P,obj.MIN_P);
            valveAmts = zeros(1,length(upstreamPVals));
            for i=1:length(upstreamPVals)
                downstreamPVals = obj.getNearbyDataPoints(downstreamPressure,min(obj.MAX_P,upstreamPVals(i)-obj.DATA_SPACING),obj.MIN_P);
%                 disp(downstreamPVals);
                innerValveAmts = zeros(1,length(downstreamPVals));
                for j=1:length(downstreamPVals)
                    innerValveAmts(j) = obj.getValvePosForMDotOfDatapoint(upstreamPVals(i),downstreamPVals(j),mdot);
                end
                if(length(downstreamPVals) < 2)
                   valveAmts(i) = innerValveAmts(1); 
                else
                    valveAmts(i) = ppval(spline(downstreamPVals,innerValveAmts),downstreamPressure);
                end
            end
            valveAmt = ppval(spline(upstreamPVals,valveAmts),upstreamPressure);
%             disp("Found valve amt: "+valveAmt);
        end
        
        function mdot = getMassFlow(obj,upstreamPressure,downstreamPressure,valveAmt)
            upstreamPVals = obj.getNearbyDataPoints(upstreamPressure,obj.MAX_P,obj.MIN_P);
            mdots = zeros(1,length(upstreamPVals));
            for i=1:length(upstreamPVals)
                downstreamPVals = obj.getNearbyDataPoints(downstreamPressure,min(obj.MAX_P,upstreamPVals(i)-obj.DATA_SPACING),obj.MIN_P);
%                 disp(downstreamPVals);
                innerMDots = zeros(1,length(downstreamPVals));
                for j=1:length(downstreamPVals)
                    innerMDots(j) = obj.getFlowRateOfDatapoint(upstreamPVals(i),downstreamPVals(j),valveAmt);
                end
                if(length(downstreamPVals) < 2)
                   mdots(i) = innerMDots(1); 
                else
                    mdots(i) = ppval(spline(downstreamPVals,innerMDots),downstreamPressure);
                end
            end
            mdot = ppval(spline(upstreamPVals,mdots),upstreamPressure);
        end
    end
    methods(Access=private)
        function valvePos = getValvePosForMDotOfDatapoint(obj,upstreamPressure,downstreamPressure,mdot)
            coeffs = obj.dataInverted([num2str(upstreamPressure),'|',num2str(downstreamPressure)]);
            coeffs(length(coeffs)) = 0; %Force zero flow for closed valves
            valvePos = polyval(coeffs,mdot);
%             disp("P1: "+upstreamPressure+", P2: "+downstreamPressure+", mdot: "+mdot+" valvePos: "+valvePos);
        end
        
        function mdot = getFlowRateOfDatapoint(obj,upstreamPressure,downstreamPressure,valveAmt)
            coeffs = obj.data([num2str(upstreamPressure),'|',num2str(downstreamPressure)]);
            coeffs(length(coeffs)) = 0; %Force zero flow for closed valves
            mdot = polyval(coeffs,valveAmt);
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
                P2 = nearest + 1.*obj.DATA_SPACING;
                P3 = nearest + 2.*obj.DATA_SPACING;
                if(P3 > maxP)
                   P3 = nearest; 
                   if(P2 > maxP)
                       P2 = nearest;
                   end
                end
            end
            if(P2 ~= nearest)
                PVals(2) = P2;
            end
            if(P3 ~= nearest)
                PVals(3) = P3;
            end
        end
    end
end