classdef BallValve < FlowCoeffFlowRestriction
    properties(SetAccess=private)
        openAmt = 1; %1 is fully open
        fullyOpenCv = 12; %Flow coefficient for when valve fully open
        CvShapePolynomialFit = [-8.20125000000565,15.4351730769403,-8.93336538463520,2.68560314686376,-0.247832167834784,0.259200466200711,0];
    end
    methods
        function obj = BallValve(fullyOpenCv,initialOpenAmt)
            obj.fullyOpenCv = fullyOpenCv;
            obj.openAmt = initialOpenAmt;
        end
        
        function setCvShapePolynomialCoeffs(obj,P)
            obj.CvShapePolynomialFit = P;
        end
        
        function setOpenAmtFraction(obj,amt)
           obj.openAmt = amt;
        end
        
        function Cv = getFlowCoefficient(obj)
            Cv = obj.fullyOpenCv .* polyval(obj.CvShapePolynomialFit,obj.openAmt);
        end
        
        function [T,P,X,v] = getDownstreamTemperaturePressureFromMassFlow(obj,mdot,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
            v = vUpstream;
            rhoUpstream = SaturatedNitrous.getDensity(XUpstream,TUpstream,PUpstream);
            rhoVUpstream = v.*rhoUpstream;
%             if(mdot == 0)
%                T = TUpstream;
%                P = PUpstream;
%                X = XUpstream;
%                return;
%             end
            TSatUpstream = SaturatedNitrous.getSaturationTemperature(PUpstream);
            liquidUpstream = fluidType == FluidType.NITROUS_LIQUID || (fluidType == FluidType.NITROUS_GENERAL && TUpstream<TSatUpstream);
            dP = obj.getPressureChangeForMassFlow(mdot,fluidType,TUpstream,PUpstream,XUpstream);
            P = PUpstream + dP;
            
            if(imag(P) ~= 0)
                disp("P upstream: "+PUpstream);
                disp("dP: "+dP);
                disp("mdot: "+mdot);
                disp("T upstream: "+TUpstream);
                disp("X Upstream: "+XUpstream);
               error('P is complex downstream of valve!');
            end
            
            if(isnan(P))
               ME = MException('BallValve:PNan','P is NaN for given mass flow (Can take any value?)');
               throw(ME);
            end
            
            if(P < 88e3)
               %disp("P TOO LOW!!"); 
               ME = MException('BallValve:DownstreamPTooLow','Downstream P too low for given mass flow');
               throw(ME);
            end
            
            %Isenthalpic
            isSaturated = false;
            if (fluidType == FluidType.NITROUS_GENERAL || fluidType == FluidType.NITROUS_LIQUID || fluidType == FluidType.NITROUS_GAS) && XUpstream ~= 1 && XUpstream ~= 0
               h1 = SaturatedNitrous.getSpecificEnthalpy(XUpstream,TUpstream,PUpstream); 
               isSaturated = true;
            else
               if(fluidType == FluidType.NITROUS_GENERAL)
                  if(XUpstream == 1)
                      h1 = FluidType.NITROUS_GAS.getSpecificEnthalpy(TUpstream,PUpstream); 
                  elseif(XUpstream == 0)
                      h1 = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(TUpstream,PUpstream); 
                  end
               else
                    h1 = fluidType.getSpecificEnthalpy(TUpstream,PUpstream); 
               end
               PVap = fluidType.getVapourPressure(TUpstream);
               if(abs((PVap - PUpstream)/PUpstream) < 0.01)
                   %Only support nitrous as a saturated fluid
                  %isSaturated = fluidType == FluidType.NITROUS_GENERAL || fluidType == FluidType.NITROUS_LIQUID || fluidType == FluidType.NITROUS_GAS; 
               end
            end
            
            if(isSaturated)
               T = SaturatedNitrous.getSaturationTemperature(P);
               X = NitrousFluidCoolProp.getProperty(FluidProperty.VAPOR_QUALITY,FluidProperty.PRESSURE,P,FluidProperty.SPECIFIC_ENTHALPY,h1);
               if(X >= 0 && X<=1)
                   rho = SaturatedNitrous.getDensity(X,T,P);
                   v = rhoVUpstream ./ rho; %Continuity
                  return; 
               end
            else
               T = fluidType.getTemperatureFromPressureEnthalpy(P,h1);
               X = XUpstream;
               rho = SaturatedNitrous.getDensity(X,T,P);
               v = rhoVUpstream ./ rho; %Continuity
            end
            PSat = SaturatedNitrous.getVapourPressure(T);
            TSat = SaturatedNitrous.getSaturationTemperature(P);
            
            %disp(" "+abs(1 - P/PSat)+" "+abs(1 - T/TSat));
            if((T>=TSat-0.0001 && liquidUpstream || T<=TSat+0.0001 && ~liquidUpstream || abs(1 - P/PSat) <= 0.0001) && (fluidType == FluidType.NITROUS_GENERAL || fluidType == FluidType.NITROUS_LIQUID || fluidType == FluidType.NITROUS_GAS))
                %Should be saturated
                X = NitrousFluidCoolProp.getProperty(FluidProperty.VAPOR_QUALITY,FluidProperty.PRESSURE,P,FluidProperty.SPECIFIC_ENTHALPY,h1);
                if(X >= 0 && X <= 1)
                    T = TSat;
                elseif (X > 1)
                    X = 1;
                    T = FluidType.NITROUS_GAS.getTemperatureFromPressureEnthalpy(P,h1);
                else
                    T = FluidType.NITROUS_LIQUID.getTemperatureFromPressureEnthalpy(P,h1);
                    X = 0;
                end
            end
            rho = SaturatedNitrous.getDensity(X,T,P);
            v = rhoVUpstream ./ rho; %Continuity
        end
        
        function [T,mdot,X,v] = getDownstreamTemperatureMassFlowFromPressureChange(obj,dP,fluidType,TUpstream,PUpstream,XUpstream,vUpstream)
            v = vUpstream;
            rhoUpstream = SaturatedNitrous.getDensity(XUpstream,TUpstream,PUpstream);
            rhoVUpstream = v.*rhoUpstream;
            TSatUpstream = SaturatedNitrous.getSaturationTemperature(PUpstream);
            liquidUpstream = fluidType == FluidType.NITROUS_LIQUID || (fluidType == FluidType.NITROUS_GENERAL && TUpstream<TSatUpstream);
            mdot = obj.getMassFlowForPressureChange(dP,fluidType,TUpstream,PUpstream,XUpstream);
            P = PUpstream+dP;
            
            %Isenthalpic
            isSaturated = false;
            if (fluidType == FluidType.NITROUS_GENERAL || fluidType == FluidType.NITROUS_LIQUID || fluidType == FluidType.NITROUS_GAS) && XUpstream ~= 1 && XUpstream ~= 0
               h1 = SaturatedNitrous.getSpecificEnthalpy(XUpstream,TUpstream,PUpstream); 
               isSaturated = true;
            else
               h1 = fluidType.getSpecificEnthalpy(TUpstream,PUpstream); 
               PVap = fluidType.getVapourPressure(TUpstream);
               if(abs((PVap - PUpstream)/PUpstream) < 0.01)
                   %Only support nitrous as a saturated fluid
                  %isSaturated = fluidType == FluidType.NITROUS_GENERAL || fluidType == FluidType.NITROUS_LIQUID || fluidType == FluidType.NITROUS_GAS; 
               end
            end
            
            if(isSaturated)
               T = SaturatedNitrous.getSaturationTemperature(P);
               X = NitrousFluidCoolProp.getProperty(FluidProperty.VAPOR_QUALITY,FluidProperty.PRESSURE,P,FluidProperty.SPECIFIC_ENTHALPY,h1);
               rho = SaturatedNitrous.getDensity(X,T,P);
               v = rhoVUpstream ./ rho; %Continuity
               if(X >= 0 && X<=1)
                  return; 
               end
            else
               T = fluidType.getTemperatureFromPressureEnthalpy(P,h1);
               X = XUpstream;
               rho = SaturatedNitrous.getDensity(X,T,P);
               v = rhoVUpstream ./ rho; %Continuity
            end
            PSat = SaturatedNitrous.getVapourPressure(T);
            TSat = SaturatedNitrous.getSaturationTemperature(P);
            
            %disp(" "+abs(1 - P/PSat)+" "+abs(1 - T/TSat));
            if((T>=TSat-0.0001 && liquidUpstream || T<=TSat+0.0001 && ~liquidUpstream || abs(1 - P/PSat) <= 0.0001) && (fluidType == FluidType.NITROUS_GENERAL || fluidType == FluidType.NITROUS_LIQUID || fluidType == FluidType.NITROUS_GAS))
                %Should be saturated
                X = NitrousFluidCoolProp.getProperty(FluidProperty.VAPOR_QUALITY,FluidProperty.PRESSURE,P,FluidProperty.SPECIFIC_ENTHALPY,h1);
                if(X >= 0 && X <= 1)
                    T = TSat;
                elseif (X > 1)
                    X = 1;
                    T = FluidType.NITROUS_GAS.getTemperatureFromPressureEnthalpy(P,h1);
                else
                    T = FluidType.NITROUS_LIQUID.getTemperatureFromPressureEnthalpy(P,h1);
                    X = 0;
                end
            end
            
            rho = SaturatedNitrous.getDensity(X,T,P);
            v = rhoVUpstream ./ rho; %Continuity
        end
    end
end