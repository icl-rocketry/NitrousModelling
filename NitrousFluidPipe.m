%Prefer to use simscape/simulink model
%
%Class to calculate required quantities associated with a flow coming from
%a reservoir out the end of a pipe driven by the pressure gradient

classdef NitrousFluidPipe < matlab.mixin.Copyable%handle %Handle class so that passed by reference not value, makes behave sanely
    properties %Allowed to be modified
        outletCrossSectionA; %Area of outlet cross section in m^3
        thinnestCrossSectionA; %Smallest cross section of pipe in m^3
    end
    properties (SetAccess=private) %Not allowed to be externally modified as may result in non-physical results
        
    end
    properties (Dependent) %Dependent properties
        
    end
    properties (Access=private)
        %Any private properties here
        upstreamConditions; %Function handle to return [Po,To] upstream
        downstreamConditions; %Function handle to return [P] downstream
    end
    properties (Constant)
       g = 9.81;
    end
    
    methods(Access=private)
        function [Po,To,fluidType] = getUpstreamTotalPressureTemperatureState(obj)
            [Po,To,fluidType] = obj.upstreamConditions();
        end
        
        function [P] = getDownstreamPressure(obj)
            P = obj.downstreamConditions(); 
        end
    end
    
    methods
        %Create a nitrous fluid pipe object.
        %Args: (upstreamConditions,downstreamConditions,outletCrossSectionA,thinnestCrossSectionA). Temp in KELVIN,
        %upstreamConditions should be a function handle that accepts no
        %params, and returns [Po,To,fluidType], where fluidType is of type
        %FluidType
        %downstreamContitions should be a function handle that accepts no
        %params, and returns [P]
        %rest in standard SI units
        function obj = NitrousFluidPipe(upstreamConditions,downstreamConditions,outletCrossSectionA,thinnestCrossSectionA) %Constructor for object
            obj.upstreamConditions = upstreamConditions;
            obj.downstreamConditions = downstreamConditions;
            obj.outletCrossSectionA = outletCrossSectionA;
            obj.thinnestCrossSectionA = thinnestCrossSectionA;
        end
        
        %Calculate the mass flow rate (Kg/s) through the pipe for current upstream
        %and downstream conditions
        function mdot = calcMassFlow(obj)
            [Po,To,fluidType] = obj.getUpstreamTotalPressureTemperatureState();
            PDownstream = obj.getDownstreamPressure();
%             disp("P Upstream:"+Po);
%             disp("P Downstream:"+PDownstream);
%             drawnow;
            if fluidType == FluidType.NITROUS_LIQUID %If liquid then treat as saturated flow
               [X2,T2,v2,h2,G] = SaturatedNitrous.getDownstreamIsentropicSaturatedHEMFlowCond(0,To,Po,PDownstream,0); 
                mdot = real(G * obj.thinnestCrossSectionA);
               return;
            end
            
            %Calculate the pressure required for choking the flow
            P = RealFlow.getEstimateOfRequiredPressureForMach(Po,To,1,fluidType);
            %disp("Approx P to choke: "+P);
            if ((PDownstream-P)/P) < 0.05 %If downstream pressure within 5% of estimate of being choked, or choked then evaluate more precisely
                try
                    [TCalc,vCalc,P,~] = RealFlow.getIsentropicConditionsForMach(Po,To,1,fluidType);
                catch
                    %Assume is not choked, as otherwise should have solved
                    P = -1;
                end
            end
            if PDownstream < P %Flow is choked
%                 disp("vCalc: "+vCalc);
                mdot = fluidType.getDensity(TCalc,P) * vCalc * obj.thinnestCrossSectionA;
%                 disp("mDot: "+mdot);
            else
                %Flow not choked, calculate the condition at the
                %end of the pipe
                CpStagnation = fluidType.getCp(To,Po);
                [TCalc,vCalc,M] = RealFlow.getIsentropicTempVelocity(Po,To,To,CpStagnation,PDownstream,fluidType);
%                 disp("vCalc: "+vCalc);
                mdot = fluidType.getDensity(TCalc,P) * vCalc * obj.outletCrossSectionA;
%                 disp("mDot: "+mdot);
            end
%             drawnow;
        end
        
        %Calculate the mass flow rate (Kg/s), temperature (K) and velocity
        %(m/s) at the outlet of this pipe for given upstream and downstream
        %conditions
        function [mdot,temp,vel] = calcOutlet(obj)
            [Po,To,fluidType] = obj.getUpstreamTotalPressureTemperatureState();
            PDownstream = obj.getDownstreamPressure();
            
            if fluidType == FluidType.NITROUS_LIQUID %If liquid then treat as saturated flow
               [X2,T2,v2,h2,G] = SaturatedNitrous.getDownstreamIsentropicSaturatedHEMFlowCond(0,To,Po,PDownstream,0); 
                mdot = G * obj.thinnestCrossSectionA;
                temp = T2;
                vel = v2;
               return;
            end
            
            %Calculate the pressure required for choking the flow
            P = RealFlow.getEstimateOfRequiredPressureForMach(Po,To,1,fluidType);
            if ((PDownstream-P)/P) < 0.1 %If downstream pressure within 10% of estimate for being choked, or is choked then evaluate more precisely
                try
                    [TCalc,vCalc,P,~] = RealFlow.getIsentropicConditionsForMach(Po,To,1,fluidType);
                catch
                    %Assume is not choked, as otherwise should have solved
                    P = -1;
                end
            end
            if PDownstream < P %Flow is choked
                chokedDensity = fluidType.getDensity(TCalc,P);
                mdot = chokedDensity * vCalc * obj.thinnestCrossSectionA;
                CpStagnation = fluidType.getCp(To,Po);
                [temp,vel] = RealFlow.getVelocityTempForKnownMassFlow(To,PDownstream,chokedDensity,vCalc,obj.thinnestCrossSectionA,obj.outletCrossSectionA,CpStagnation,fluidType);
            else
                %Flow not choked, calculate the condition at the
                %end of the pipe
                CpStagnation = fluidType.getCp(To,Po);
                [TCalc,vCalc,M] = RealFlow.getIsentropicTempVelocity(Po,To,To,CpStagnation,PDownstream,fluidType);
                %fprintf(['CpStag: ',num2str(CpStagnation),' vCalc: ',num2str(vCalc),'\n']);
                mdot = fluidType.getDensity(TCalc,P) * vCalc * obj.outletCrossSectionA;
                temp = TCalc;
                vel = vCalc;
            end
        end
    end
end