% NOTE! If TankEnergyBalanceFast.USE is true then this will lookup based
% on the precomputed values! (Faster)
%
% Script to model the state of a nitrous tank with constant cross sectional area, with given pertubations to mass etc.
% NOTE this model assumes all fluids within tank are always at the same
% temperature and in equilibrium so will transfer heat between them at infinite rate so as to
% keep it always this way.
% NOTE this model is simplified and essentially pretends air does not
% exist. The only correction related to this made is that it arbitrarily
% enforces nitrous vapour's pressure must be always >= atmospheric
% pressure. (Eg. so that nitrous density doesn't have it near vacuum for
% low nitrous mass)
% Depends on SaturatedNitrous.m,NonLinearEquationSolver.m,NitrousFluid.m
% By Eddie Brown

classdef GeometricNitrousTank < matlab.mixin.Copyable%handle %Handle class so that passed by reference not value, makes behave sanely
    properties %Allowed to be modified
        
    end
    properties (SetAccess=private) %Not allowed to be externally modified as may result in non-physical results
        mTotalNitrous; %Total mass of nitrous contents of the tank  
        tankTotalVolume; %Tank's total volume
        crossSectionArea; %Cross sectional area of the tank (m^2)
        temp; %Current temperature of the contents of the tank, IN KELVIN. Single temp for all contents since heat transfer in saturated state between liquid and vapour very large and as flow into and out of tank will induce large amounts of convection
        gasVentHoleHeight;
    end
    properties (Dependent) %Dependent properties
        mTotalTank; %Current total mass in tank
        tankHeight; %Height of the tank (in metres)
        mLiquid; %Current mass of liquid in the tank
        vapourPressure; %The current vapour pressure of the vapour in the tank
        mVapour; %Current mass of nitrous vapour in the tank
        volVapour; %Current volume of the vapour in the tank
        volLiquid; %Current volume of the liquid in the tank
        liquidDensity; %Current density of the liquid in the tank
        vapourDensity; %Current density of the vapour in the tank
        liquidHeight; %Current height of the liquid nitrous (m)
        vapourHeight; %Current height of the nitrous vapour (m)
        pressureAtBaseOfTank; %Current pressure at the base of the tank
        isSaturated; %True if nitrous is in saturated condition
    end
    properties (Access=private)
        %Any private properties here
    end
    properties (Constant)
       %Constants here 
       %The amount of liquid for which drainSmallAmountOfLiquid is
       %considered accurate. Smaller will be more accurate but slower to
       %run calculations of mass drainage effect, larger will be less
       %accurate but faster to run calculations of mass drainage effect
       PREFERRED_SMALL_DRAIN_AMOUNT = 0.001; %1 grams
       R_AIR = 287;
       Cp_AIR = 1005;
       g = 9.81;
       ATMOSPHERIC_PRESSURE = 101325;
    end
    
    methods(Access=private)        
        
        %Define getter function for mass of liquid in tank
        function mLiq=getSaturatedLiquidMass(obj)
            vDensity = SaturatedNitrous.getVapourDensity(obj.temp);
            lDensity = SaturatedNitrous.getLiquidDensity(obj.temp);
            %From trivial re-arrangement of volume and mass summations
            mLiq = (obj.tankTotalVolume - (obj.mTotalNitrous/vDensity)) / ((1/lDensity)-(1/vDensity));
        end
        
        %Define getter function for mass of vapour in tank
        function mVap=getSaturatedVapourMass(obj)
            vDensity = SaturatedNitrous.getVapourDensity(obj.temp);
            lDensity = SaturatedNitrous.getLiquidDensity(obj.temp);
            %From trivial re-arrangement of volume and mass summations
            mVap = (obj.tankTotalVolume - (obj.mTotalNitrous/lDensity)) / ((1/vDensity)-(1/lDensity));
        end
    end
    
    methods
        %Create a nitrous tank object.
        %Args: (temp, mNitrousTotal, mAir,tankHeight,crossSectionalArea). Temp in KELVIN,
        %rest in standard SI units
        function obj = GeometricNitrousTank(temp,mNitrousTotal,tankHeight,crossSectionalArea,gasVentHoleHeight) %Constructor for object
            obj.tankTotalVolume = tankHeight*crossSectionalArea; %Set volume of tank
            obj.mTotalNitrous = mNitrousTotal;
            obj.crossSectionArea = crossSectionalArea;
            obj.temp = temp; %Set temp of tank
            obj.gasVentHoleHeight = gasVentHoleHeight;
        end
        
        function E = getInternalEnergy(obj)
%             E1 = obj.mAir * FluidType.AIR.getSpecificInternalEnergy(obj.temp,obj.airPressure);
%             if obj.mAir == 0
%                E1 = 0; %Stops equalling NaN
%             end
            E1 = 0;
            %PLiq = obj.getPressureAtHeight(0.5*obj.liquidHeight);
            if(abs(obj.mLiquid) < 0.0001)
                E2 = 0;
            else
                uL = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_INTERNAL_ENERGY,FluidProperty.TEMPERATURE,obj.temp,FluidProperty.PRESSURE,obj.vapourPressure);
                E2 = obj.mLiquid * uL; %FluidType.NITROUS_LIQUID.getSpecificInternalEnergy(obj.temp,obj.vapourPressure);
            end
            if(abs(obj.mVapour) < 0.0001)
                E3 = 0;
            else
                uG = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.SPECIFIC_INTERNAL_ENERGY,FluidProperty.TEMPERATURE,obj.temp,FluidProperty.PRESSURE,obj.vapourPressure);
                E3 = obj.mVapour * uG;%FluidType.NITROUS_GAS.getSpecificInternalEnergy(obj.temp,obj.vapourPressure);
            end
            E = E1 + E2 + E3;
        end
        
        function E = getTankTotalEnthalpy(obj)
%             E1 = obj.mAir * FluidType.AIR.getSpecificEnthalpy(obj.temp,obj.airPressure);
%             if obj.mAir == 0
%                E1 = 0; %Stops equalling NaN
%             end
            E1 = 0;
            %PLiq = obj.getPressureAtHeight(0.5*obj.liquidHeight);
            if(obj.mLiquid == 0)
                E2 = 0; 
            else
                E2 = obj.mLiquid * FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.temp,obj.vapourPressure);
            end
            E3 = obj.mVapour * FluidType.NITROUS_GAS.getSpecificEnthalpy(obj.temp,obj.vapourPressure);
            E = E1 + E2 + E3;
        end
        
        %Sets the temperature of the tank to a given value ignoring the heat transfer required to do it, and then
        %updates the nitrous/air ratio of the tank (Which may alter the
        %temperature to be slightly different)
        function forceSetTemp(obj,temp)
            obj.temp = temp;
        end
        
        %Heat the tank by a certain amount (J), use negative heat for
        %cooling
        function addHeat(obj,heatAmount)
            before = copy(obj);

            %Internal energy in the tank before heat addition
            ECv1 = before.getInternalEnergy();
            dE = heatAmount; %Energy added (The heat)
            ECv2 = ECv1 + dE; %New internal energy of tank
            
            %Solve for temp that gives final mVapour and mLiquid with
            %correct internal energy
            solution = betterfzero(@errFunc,obj.temp,200,SaturatedNitrous.T_CRIT-1,2);
            T2 = real(solution(1));
            obj.temp = T2;
            
            function err = errFunc(endTemp)
                endTemp = real(endTemp);
                obj.temp = endTemp; %Changes masses of vapour and liquid within CV
                
                E2 = obj.getInternalEnergy();
                err = E2 - ECv2;
            end
        end
        
        function forceSetNitrousMass(obj,mNitrous)
            obj.mTotalNitrous = mNitrous;
        end
        
        %Calculate the rate of change of total internal energy (not specific) of the tank
        %required for the tank to stay at constant temperature for a mass
        %only pertubation of netMdotFill per second (Eg. if after 1 second
        %the total mass of nitrous in the tank has increased by
        %netMdotFill Kg, what will the tank's internal energy have to have
        %changed by to still have the tank at the same temp)
        function dEdt = findIntEnergyChangeRateForConstTemperatureWithFillRate(obj,netMdotFill)
            clone = copy(obj);
            dt = 1e-3; %1 millis
            clone.mTotalNitrous = obj.mTotalNitrous + netMdotFill.*dt; %Perturb mass only
            dE = clone.getInternalEnergy() - obj.getInternalEnergy(); %Change in tank internal energy
            dEdt = dE./dt;
        end
        
        function dEdt = findIntEnergyChangeRateForTempChangeRateWithFillRate(obj,dTdt,netMdotFill)
            clone = copy(obj);
            dt = 1e-2;
            clone.mTotalNitrous = obj.mTotalNitrous + netMdotFill.*dt;
            clone.temp = clone.temp + dTdt.*dt;
            dE = clone.getInternalEnergy() - obj.getInternalEnergy();
            dEdt = dE/dt;
        end
        
        function Q = findHeatRateInputForTempChangeRateWithLiquidDrainingRate(obj,dTdt,mdotLiq)
            try
                h = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.temp,obj.vapourPressure);
            catch
                h = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,obj.temp,FluidProperty.VAPOR_QUALITY,0);
            end
            H = h*mdotLiq; %Energy being lost due to liquid
            dEdt = obj.findIntEnergyChangeRateForTempChangeRateWithFillRate(dTdt,-mdotLiq);
            Q = dEdt + H;
        end
        
        function Q = findHeatRateInputToKeepTempForDrainingLiquidRate(obj,mdotLiq)
            h = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.temp,obj.vapourPressure);
            H = h*mdotLiq; %Energy being lost due to liquid
            dEdt = obj.findIntEnergyChangeRateForConstTemperatureWithFillRate(-mdotLiq); %dEdt to keep temp same
            Q = dEdt + H;
        end
        
        function Q = findHeatToMakeTemp(obj,tempDesired)
            clone = copy(obj);
            clone.temp = tempDesired;
            Q = clone.getInternalEnergy() - obj.getInternalEnergy();
        end
        
        %Change the nitrous in the tank's mass and energy. Positive for
        %energy/mass input, and negative for energy/mass output. EIn is not
        %specific, but should take into account mass and type of each mass
        %in it's calculation
        function setNitrousMassEnergyAndAdjustTemp(obj,m,E)
            obj.mTotalNitrous = m;
            %Solving for in saturated state max T until modelling can do
            %supercritical
            T2 = betterfzero(@errTemp,obj.temp,182.23,SaturatedNitrous.T_CRIT-2,2); %Solve for T that gives correct internal energy
            %Update tank temp and nitrous vol
            obj.temp = real(T2);
            
            function err = errTemp(T)       
               obj.temp = real(T);
               if(obj.temp < 182.23) %Temp too low
                   err = -1000000;
                   return;
               end
               if(obj.temp > SaturatedNitrous.T_CRIT)
                   err = 1000000;
                   return;
               end
               E2 = obj.getInternalEnergy();               
               err = E2 - E;
%                fprintf(['T: ',num2str(T),', E: ',num2str(E2), ' Err: ',num2str(err),'\n']);
            end
        end
        
        %Change the nitrous in the tank's mass and energy. Positive for
        %energy/mass input, and negative for energy/mass output. EIn is not
        %specific, but should take into account mass and type of each mass
        %in it's calculation
        function changeNitrousMassEnergy(obj,mIn,EIn)
            %Internal energy in the tank before
            ECv1 = obj.getInternalEnergy();
            %disp("E Before: "+ECv1);
%             disp("E change: "+EIn);
%             fprintf(['E1: ',num2str(ECv1),', E2: ',num2str(ECv2),'\n']);
%             fprintf(['M1: ',num2str(obj.mTotalNitrous),', M2: ',num2str(obj.mTotalNitrous + mIn),'\n']);
            obj.mTotalNitrous = obj.mTotalNitrous + mIn; %Nitrous mass of tank after
            
            prevTemp = obj.temp;
            nitrousDensity = obj.mTotalNitrous / obj.tankTotalVolume;
            specificE = (ECv1 + EIn) / obj.mTotalNitrous;
            if(TankEnergyBalanceFast.USE)
                obj.temp = TankEnergyBalanceFast.getInstance().getTemp(nitrousDensity,specificE);
                if(obj.temp < 182.23)
                   obj.temp = 182.23; 
                end
                if(~isnan(obj.temp))
                    return;
                end
                %disp("NOT found temp for nitrous density: "+nitrousDensity+" with specific E "+(specificE)+" nitrous mass: "+obj.mTotalNitrous);
            end
            if(obj.mTotalNitrous <= 0)
                obj.temp = SaturatedNitrous.getSaturationTemperature(101325);
                return;
            end
            ECv2 = ECv1 + EIn; %Internal energy of the tank after
            clone = copy(obj);
            clone.mTotalNitrous = obj.mTotalNitrous;
            %Solving for in saturated state max T until modelling can do
            %supercritical
            try
                if prevTemp <= 184.7 &&  specificE <= 339000
                    T2 = 182.23;
                else
                    T2 = betterfzero(@errTemp,prevTemp,182.23,SaturatedNitrous.T_CRIT-2,2,1000,0.5e-50,true); %Solve for T that gives correct internal energy
                end
            catch err
                %Initial tank conditions may disagree slightly with energy
%                 disp(prevTemp+" "+nitrousDensity+" "+specificE);
                %nitrousDensity <= 0.03 &&
                if prevTemp <= 184.7 &&  specificE <= 339000
                    T2 = 182.23;
                else
                    rethrow(err);
                end
            end
            
            %Update tank temp and nitrous vol
            obj.temp = real(T2);
            
            function err = errTemp(T)       
               clone.temp = real(T);
               if(clone.temp < 182.23) %Temp too low
                   err = -1000000;
                   return;
               end
               if(clone.temp > SaturatedNitrous.T_CRIT)
                   err = 1000000;
                   return;
               end
               E2 = clone.getInternalEnergy();               
               err = E2 - ECv2;
%                fprintf(['T: ',num2str(T),', E: ',num2str(E2), ' Err: ',num2str(err),'\n']);
            end
        end
        
        %Adds small amount of nitrous to the tank. EIn should be all energy input as a result of this mass addition, it should NOT be specific (NOT per Kg) 
        function addAmountOfNitrous(obj,mIn,EIn)
            obj.changeNitrousMassEnergy(mIn,EIn);
        end
        
        %Adds small amount of nitrous to the tank. EIn should be all energy input as a result of this mass addition, it should NOT be specific (NOT per Kg) 
        function addSmallAmountOfNitrous(obj,mIn,EIn)
            obj.changeNitrousMassEnergy(mIn,EIn);
        end
        
        %Re-calculates the state of the tank after draining an amount of
        %liquid from it.
        %Works by splitting amount into tiny amounts and then calling
        %drainSmallAmountOfLiquid(mDrained) until correct amount drained
        function drainAmountOfLiquid(obj,mDrained)
            h = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.temp,obj.getPressureAtHeight(0));
            %h = FluidType.NITROUS_LIQUID.getSpecificInternalEnergy(obj.temp,obj.getPressureAtHeight(0));
            H = h*mDrained;
            obj.changeNitrousMassEnergy(-mDrained,-H);
        end
        
        function drainSmallAmountOfLiquid(obj,mDrained)
           h = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.temp,obj.getPressureAtHeight(0));
           %h = FluidType.NITROUS_LIQUID.getSpecificInternalEnergy(obj.temp,obj.getPressureAtHeight(0));
           H = h*mDrained;
           obj.changeNitrousMassEnergy(-mDrained,-H);
        end
        
        %Re-calculates the state of the tank after draining an amount of
        %vapour from it.
        %Works by splitting amount into tiny amounts and then calling
        %drainSmallAmountOfVapour(mDrained) until correct amount drained
        function drainAmountOfVapour(obj,mDrained)
            h = FluidType.NITROUS_GAS.getSpecificEnthalpy(obj.temp,obj.getPressureAtHeight(obj.gasVentHoleHeight));
            H = h*mDrained;
            obj.changeNitrousMassEnergy(-mDrained,-H);
        end
        
        function h = specificEnthalpyVapour(obj)
            h = FluidType.NITROUS_GAS.getSpecificEnthalpy(obj.temp,obj.getPressureAtHeight(obj.gasVentHoleHeight));
        end
        
        function drainSmallAmountOfVapour(obj,mDrained)
           h = FluidType.NITROUS_GAS.getSpecificEnthalpy(obj.temp,obj.getPressureAtHeight(obj.gasVentHoleHeight));
           H = h*mDrained;
           obj.changeNitrousMassEnergy(-mDrained,-H);
        end
        
        function P = getPressureAtHeight(obj,h)
           P = obj.vapourPressure;
           if(~obj.isSaturated())
              %Assume calculated P is at half tank height
%              P = P - obj.vapourDensity * GeometricNitrousTank.g * (h-(obj.tankHeight/2));
              return;
           end
           liqHeight = obj.liquidHeight; %Height of nitrous liquid phase
           if h < liqHeight
              P = P + obj.liquidDensity * GeometricNitrousTank.g * (liqHeight-h);
           else
              P = P - obj.vapourDensity * GeometricNitrousTank.g * (h-liqHeight);
           end
        end
        
        function mTotalTank = get.mTotalTank(obj)
           mTotalTank = obj.mTotalNitrous; 
        end
        
        %Define getter function for mass of liquid in tank
        function mLiq=get.mLiquid(obj)
            if ~obj.isSaturated() %If not at saturation conditions, then assume all nitrous is a gas
               if(obj.vapourPressure > SaturatedNitrous.getVapourPressure(obj.temp))
                    mLiq = obj.mTotalNitrous;
               else
                   mLiq = 0;
               end
               return;
            end
            %From trivial re-arrangement of volume and mass summations
            mLiq = obj.getSaturatedLiquidMass();
%             if (mLiq < 0)
%                warning("Liquid mass negative, no physical saturated solution can exist for these conditions"); 
%             end
        end
        
        %Define getter function for mass of vapour in tank
        function mVap=get.mVapour(obj)
            if ~obj.isSaturated() %If not at saturation conditions, then assume all nitrous is a gas
                if(obj.vapourPressure < SaturatedNitrous.getVapourPressure(obj.temp))
                    mVap = obj.mTotalNitrous;
                else
                   mVap = 0; 
                end
               return;
            end
            %From trivial re-arrangement of volume and mass summations
            mVap = obj.getSaturatedVapourMass();
            if (mVap < 0)
               warning("Vapour mass negative, no physical saturated solution can exist for these conditions"); 
            end
        end
        
        %Define getter function for height of the tank
        function val=get.tankHeight(obj)
           val = obj.tankTotalVolume / obj.crossSectionArea; 
        end
        
        %Define getter function for current pressure exerted by vapour in
        %the tank
        function val=get.vapourPressure(obj)
            if ~obj.isSaturated()
                if obj.mTotalNitrous == 0 %If no nitrous in the tank
                   val = obj.ATMOSPHERIC_PRESSURE; %Just return air pressure
                   return;
                end
                
                %Calculate pressure for if tank was entirely filled by
                %nitrous
               rho = obj.mTotalNitrous / obj.tankTotalVolume;
%                disp("R: "+rho+" "+obj.mTotalNitrous+" "+obj.tankTotalVolume);
%                drawnow;
               val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.PRESSURE,FluidProperty.TEMPERATURE,obj.temp,FluidProperty.DENSITY,rho);
%                disp("P tank: "+val+" sat m liq: "+obj.getSaturatedLiquidMass()+" sat m vap: "+obj.getSaturatedVapourMass());
               if val > SaturatedNitrous.getVapourPressure(obj.temp)
                   val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.PRESSURE,FluidProperty.TEMPERATURE,obj.temp,FluidProperty.DENSITY,rho);
               end
               %If mass of nitrous is too low and pressure is not physical
               %(less than atmospheric)
               if(val < obj.ATMOSPHERIC_PRESSURE)
                    val = obj.ATMOSPHERIC_PRESSURE;
               end
               return;
            end
            val = SaturatedNitrous.getVapourPressure(obj.temp);
        end
        
        %Define getter function for current liquid density
        function val=get.liquidDensity(obj)
            val = SaturatedNitrous.getLiquidDensity(obj.temp);
        end
        
        %Define getter function for current vapour density
        function val=get.vapourDensity(obj)
            if ~obj.isSaturated() %If not saturated, then get density of nitrous gas
                %is just m/V
%                 disp("Temp: "+obj.temp);
%                 disp("Taking pressure to be "+obj.vapourPressure);
               val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.DENSITY,FluidProperty.TEMPERATURE,obj.temp,FluidProperty.PRESSURE,obj.vapourPressure);
               return;
            end
            val = SaturatedNitrous.getVapourDensity(obj.temp);
        end
        
        %Define getter function for current liquid volume
        function vol=get.volLiquid(obj)
            if ~obj.isSaturated() %If not saturated assume all nitrous is vapour
               vol = 0; 
               return;
            end
           vol = obj.mLiquid / obj.liquidDensity; %V = m/rho
        end
        
        %Define getter function for current vapour volume
        function vol=get.volVapour(obj)
%             if ~obj.isSaturated() %If not saturated then vapour takes up entire tank
%                vol = obj.mTotalNitrous ./ obj.vapourDensity;
%                return;
%             end
           vol = obj.mVapour / obj.vapourDensity; %V = m/rho
        end
        
        %Define getter function for current liquid height
        function val=get.liquidHeight(obj)
            val = obj.volLiquid / obj.crossSectionArea;
        end
        
        %Define getter function for current vapour height
        function val=get.vapourHeight(obj)
            val = obj.volVapour / obj.crossSectionArea;
        end
        
        function val = get.isSaturated(obj)
            %If both >=0 then must be saturated, if one is negative then
            %not - instead assume nitrous is all vapour
%             if true
%                 %TODO Fix without breaking vapour pressure calculations
%                 val = true;
%                 return;
%             end
           val = obj.getSaturatedLiquidMass() >= 0 && obj.getSaturatedVapourMass() >= 0;
        end
        
        %Define getter function for the pressure at the base of the tank
        function val=get.pressureAtBaseOfTank(obj)
            %P1 = P + rho*g*h
%             if obj.mTotalNitrous == 0
%                 val = obj.airPressure + obj.airDensity * GeometricNitrousTank.g * (obj.airHeight/2);
%                return; 
%             end
            if(obj.isSaturated())
                val = obj.vapourPressure + obj.liquidDensity * GeometricNitrousTank.g * obj.liquidHeight;
            end
            %TODO
        end
    end
end
