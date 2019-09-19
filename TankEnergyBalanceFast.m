classdef TankEnergyBalanceFast < matlab.mixin.Copyable
    properties
        data; %The pre baked data
    end
    properties  (Constant)
        TANK_VOL_USED_FOR_CALCS = 0.01767;
        DATA_SPACING_E = 500;
        MAX_E = 4.2e5;
        MIN_E = 0;
        DATA_SPACING_NITROUS_DENSITY = 0.010 / TankEnergyBalanceFast.TANK_VOL_USED_FOR_CALCS;
        MIN_NITROUS_DENSITY = 0;
        MAX_NITROUS_DENSITY = 20 / TankEnergyBalanceFast.TANK_VOL_USED_FOR_CALCS;
        USE = true;
    end
    methods(Static)
        function bal = getInstance()
           persistent instance;
%             instance = [];
           if(isempty(instance))
              instance = TankEnergyBalanceFast('preBakedData/nitrousTankInternalEnergies.mat');
           end
           bal = instance;
           return;
        end
    end
    methods
        function obj = TankEnergyBalanceFast(datafileName)
            obj.data = load(datafileName,'-mat');
            obj.data = obj.data.dataMap; %Extract map from struct
        end
        
%         function mdot = getMassFlow(obj,upstreamPressure,downstreamPressure,valveAmt)
%             upstreamPVals = obj.getNearbyDataPoints(upstreamPressure,obj.MAX_P,obj.MIN_P);
%             mdots = zeros(1,length(upstreamPVals));
%             for i=1:length(upstreamPVals)
%                 downstreamPVals = obj.getNearbyDataPoints(downstreamPressure,min(obj.MAX_P,upstreamPVals(i)-obj.DATA_SPACING),obj.MIN_P);
%                 %                 disp(downstreamPVals);
%                 innerMDots = zeros(1,length(downstreamPVals));
%                 for j=1:length(downstreamPVals)
%                     innerMDots(j) = obj.getFlowRateOfDatapoint(upstreamPVals(i),downstreamPVals(j),valveAmt);
%                 end
%                 if(length(downstreamPVals) < 2)
%                     mdots(i) = innerMDots(1);
%                 else
%                     mdots(i) = ppval(spline(downstreamPVals,innerMDots),downstreamPressure);
%                 end
%             end
%             mdot = ppval(spline(upstreamPVals,mdots),upstreamPressure);
%         end
    end
    methods
        function T = getTemp(obj,massPerVol,specificE)
            nearest = round(massPerVol/obj.DATA_SPACING_NITROUS_DENSITY) * obj.DATA_SPACING_NITROUS_DENSITY;
            if(nearest > obj.MAX_NITROUS_DENSITY)
               nearest = obj.MAX_NITROUS_DENSITY; 
            end
            if(nearest < obj.MIN_NITROUS_DENSITY)
               nearest = obj.MIN_NITROUS_DENSITY;
            end
            T1 = obj.getTempForMassDatapointAndGivenE(nearest,specificE);
            if isnan(T1) %Datapoint not in map
                nearest = nearest + obj.DATA_SPACING_NITROUS_DENSITY;
                T1 = obj.getTempForMassDatapointAndGivenE(nearest,specificE);
            end
            
            M2 = nearest + obj.DATA_SPACING_NITROUS_DENSITY;
            M3 = nearest - obj.DATA_SPACING_NITROUS_DENSITY;
            TsAssigned = false;
            T2 = obj.getTempForMassDatapointAndGivenE(M2,specificE);
            if isnan(T2) %No data point to the right of T1
                M2 = M3;
                M3 = nearest - 2.*obj.DATA_SPACING_NITROUS_DENSITY;
                T2 = obj.getTempForMassDatapointAndGivenE(M2,specificE);
                T3 = obj.getTempForMassDatapointAndGivenE(M3,specificE);
%                 isnan(T3)
%                     T3 = NaN; %Not enough for three datapoints
%                 end
                TsAssigned = true;
            end
            if ~TsAssigned
                T3 = obj.getTempForMassDatapointAndGivenE(M3,specificE);
                if isnan(T3) %NO data point to the left of T1
                    M3 = nearest + 2.*obj.DATA_SPACING_NITROUS_DENSITY;
                    T3 = obj.getTempForMassDatapointAndGivenE(M3,specificE);
%                     if isnan(T3)
%                         T3 = NaN; %Not enough for three datapoints
%                     end
                end
            end
            
            if isnan(T1) || isnan(T2)
               T = NaN;
               return;
            end
            Ts = [T1,T2];
            Ms = [nearest,M2];
            if(~isnan(T3))
               Ms = [Ms,M3];
               Ts = [Ts,T3]; 
            end
            
            T = spline(Ms,Ts,massPerVol);
        end
        
        function T = getTempForMassDatapointAndGivenE(obj,massDataPoint,E)
            nearest = round(E/obj.DATA_SPACING_E) * obj.DATA_SPACING_E;
            if(nearest > obj.MAX_E)
                nearest = obj.MAX_E;
            end
            if(nearest < obj.MIN_E)
                nearest = obj.MIN_E;
            end
            T1 = obj.getTempOfDataPoint(nearest,massDataPoint);  
            if isnan(T1) %Datapoint not in map
                nearest = nearest + obj.DATA_SPACING_E;
                T1 = obj.getTempOfDataPoint(nearest,massDataPoint);
            end
            E2 = nearest + obj.DATA_SPACING_E;
            E3 = nearest - obj.DATA_SPACING_E;
            TsAssigned = false;
            T2 = obj.getTempOfDataPoint(E2,massDataPoint);
            if isnan(T2) %No data point to the right of T1
                E2 = E3;
                E3 = nearest - 2.*obj.DATA_SPACING_E;
                T2 = obj.getTempOfDataPoint(E2,massDataPoint);
                T3 = obj.getTempOfDataPoint(E3,massDataPoint);
%                 if(isnan(T3))
%                     T3 = NaN; %Not enough for three datapoints
%                 end
                TsAssigned = true;
            end
            if ~TsAssigned
                T3 = obj.getTempOfDataPoint(E3,massDataPoint);
                if isnan(T3) %NO data point to the left of T1
                    E3 = nearest + 2.*obj.DATA_SPACING_E;
                    T3 = obj.getTempOfDataPoint(E3,massDataPoint);
%                     if isnan(T3)
%                         T3 = NaN; %Not enough for three datapoints
%                     end
                end
            end
            if isnan(T1) || isnan(T2)
               T = NaN;
               return;
            end
            Ts = [T1,T2];
            Es = [nearest,E2];
            if(~isnan(T3))
               Es = [Es,E3];
               Ts = [Ts,T3]; 
            end
            
            T = spline(Es,Ts,E);
        end
        
        function T = getTempOfDataPoint(obj,tankSpecificInternalEnergy,mNitrous)
            k = [num2str(tankSpecificInternalEnergy),'|',num2str(mNitrous)];
            if ~isKey(obj.data,k)
               T = NaN;
               return;
            end
            T = obj.data(k);
        end
    end
end