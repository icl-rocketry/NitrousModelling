% Script to model the state of a nitrous tank with constant cross sectional area, with given pertubations to mass etc.
% NOTE this model assumes all fluids within tank are always at the same
% temperature and in equilibrium so will transfer heat between them at infinite rate so as to
% keep it always this way.
% Depends on SaturatedNitrous.m,NonLinearEquationSolver.m,NitrousFluid.m
% NOTE THE MODELLING OF AIR WITH THIS MODEL IS BROKEN, DO NOT USE WITH
% NON-ZERO AIR MASS
% By Eddie Brown

classdef GeometricNitrousAndAirTank < matlab.mixin.Copyable%handle %Handle class so that passed by reference not value, makes behave sanely
    properties %Allowed to be modified
        
    end
    properties (SetAccess=private) %Not allowed to be externally modified as may result in non-physical results
        mAir; %Current mass of air in the tank
        mTotalNitrous; %Total mass of nitrous contents of the tank  
        tankTotalVolume; %Tank's total volume
        volNitrous; %Curent volume of the tank taken up by nitrous - Technically a dependent property, but requires solving simultanous equations numerically to find so will make a private variable instead
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
        volAir; %Current volume of the air in the tank
        liquidDensity; %Current density of the liquid in the tank
        vapourDensity; %Current density of the vapour in the tank
        airDensity; %Current density of the air in the tank
        airPressure; %Current pressure of the air in the tank
        liquidHeight; %Current height of the liquid nitrous (m)
        vapourHeight; %Current height of the nitrous vapour (m)
        airHeight; %Current height of the air (m)
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
    end
    
    methods(Access=private)
        %Change the volume of nitrous in the tank in a thermodynamically
        %consistent way
        function setVolNitrous(obj,volNitrous)
            
           if ~obj.isSaturated() %If not saturated then nitrous is all vapour so vapour mass won't be changing so can do this in this model where temp is same for all fluids in container
              obj.volNitrous = volNitrous;
              return;
           end
           TGuess = obj.temp;
           %Internal energy of whole tank before (Should equal internal
           %energy after) - as no heat input, no external boundary work, no shaft work
           E1 = obj.getInternalEnergy();
           %Solve for temp where energy is conserved
           solution = betterfzero(@errFunc,TGuess,220,378,1*10^-7,1000);
           T2 = real(solution(1));
           obj.volNitrous = volNitrous;
           obj.temp = T2;
           
            function err = errFunc(TGuess)
                clone = copy(obj);
                clone.volNitrous = volNitrous;
                clone.temp = real(TGuess);
                
                T2 = real(TGuess);
                E2 = clone.getInternalEnergy();
                %fprintf(['Vol Nit: ',num2str(volNitrous),' T: ',num2str(clone.temp),' Err in E: ',num2str((E2-E1)),'\n']);
                err = E2-E1;
            end
        end
        
%         %Inner loop of the "drainSmallAmountOfLiquid"'s iterative solver
%         %Uses guess for mass vapourized to calculate end state final
%         %temperature and updates value for endState's temperature
%         %accordingly.
%         %This is then used to find the final liquid mass which is then used
%         %to find the amount of liquid that must have vapourized
%         %The error between the vapour amount guess and that calculated is
%         %returned so that the iterative solver can find the solution
%         %Is moved here as a separate function so can be called twice each
%         %loop without code duplication (so gradient can be approximated for
%         %newton raphson solver)
%         function error = calcVapourAmountGuess(obj,mDrained,endState,mVap)
%             mVap = real(mVap);
%             mDrained = real(mDrained);
%             CvAir = GeometricNitrousAndAirTank.Cp_AIR - GeometricNitrousAndAirTank.R_AIR;
%             %CpAir = GeometricNitrousAndAirTank.Cp_AIR;
%             if obj.isSaturated()
%                 CpVap1 = SaturatedNitrous.getVapourIsobaricHeatCapacity(obj.temp);
%                 CpLiq1 = SaturatedNitrous.getLiquidIsobaricHeatCapacity(obj.temp);
%             else 
%                 CpVap1 = NitrousFluid.getGasCp(obj.temp,obj.vapourPressure);
%                 CpLiq1 = NitrousFluid.getLiquidCp(obj.temp,obj.vapourPressure);
%             end
%             gammaVap1 = NitrousFluid.getGasSpecificHeatRatio(obj.temp,obj.vapourPressure);
%             gammaLiq1 = NitrousFluid.getLiquidSpecificHeatRatio(obj.temp,obj.vapourPressure);
%             CvVap1 = CpVap1 / gammaVap1;
%             CvLiq1 = CpLiq1 / gammaLiq1;
%             %Commented stuff left for help making future improvement to
%             %unlinking temps of tank contents
%             
%             %Calculate dQ, the heat removed from liquid during this vapourization
% %             dQ = mVap * SaturatedNitrous.getLiquidEnthalpyOfVaporisation(obj.temp);
%             
%             %First let's remove the vapour mass that left the tank
%             dE = (-mDrained) * CpLiq1 * obj.temp; %-Energy lost from vapour due to mass exiting tank (=mCpTo)
%             mVNew = obj.mVapour + mVap; %New mass of vapour now some lost
%             
%             %Now let's calculate what happened to the liquid to vapourize
%             %mVap (in response to gas leaving)
%             %Calculate dT, the temp drop of the remaining liquid nitrous
%             %due to losing this heat (From vapourising). (dT = dQ/(m*C))
% %             dT = -dQ / ((obj.mLiquid) * CpLiq1);
% %             T2 = obj.temp +dT; %Calculate end temp of liquid
%             
%             mLNew = obj.mLiquid - mVap - mDrained; %New mass of liquid
%             %Energy of contents of tank            
%             E1 = obj.mAir*CvAir*obj.temp + obj.mLiquid*CvLiq1*obj.temp + obj.mVapour*CvVap1*obj.temp + dE;
%             
%             %Uniform temperature of tank contents that gives same internal
%             %energy
%             T2 = E1 / (obj.mAir*CvAir + mLNew*CvLiq1 + mVNew*CvVap1);
%             
%             endState.temp = T2; %Update end temp of the end state
%             %Calculate mass vapourized in reaching end state, before liquid mass
%             %minus final liquid mass. Should equal mVap if mVap was correct
%             mVapourizedToReachEndState = (obj.mLiquid) - endState.mLiquid;
%             error = -(mVapourizedToReachEndState-mVap); %If positive then mVap was too big
%             %fprintf(['Vap: ',num2str(mVap),', err: ',num2str(error),' req Vap: ',num2str(mVapourizedToReachEndState),', T: ',num2str(T2),'\n']);
%         end
        
%         %Inner loop of the "drainSmallAmountOfValpour"'s iterative solver
%         %Uses guess for mass vapourized to calculate end state final
%         %temperature and updates value for endState's temperature
%         %accordingly.
%         %This is then used to find the final liquid mass which is then used
%         %to find the amount of liquid that must have vapourized
%         %The error between the vapour amount guess and that calculated is
%         %returned so that the iterative solver can find the solution
%         %Is moved here as a separate function so can be called twice each
%         %loop without code duplication (so gradient can be approximated for
%         %newton raphson solver)
%         function error = calcVapourAmountGuessForDrainVapour(obj,mDrained,endState,mVap)
%             mVap = real(mVap);
%             mDrained = real(mDrained);
%             CvAir = GeometricNitrousAndAirTank.Cp_AIR - GeometricNitrousAndAirTank.R_AIR;
%             %CpAir = GeometricNitrousAndAirTank.Cp_AIR;
%             if obj.isSaturated()
%                 CpVap1 = SaturatedNitrous.getVapourIsobaricHeatCapacity(obj.temp);
%                 CpLiq1 = SaturatedNitrous.getLiquidIsobaricHeatCapacity(obj.temp);
%             else 
%                 CpVap1 = NitrousFluid.getGasCp(obj.temp,obj.vapourPressure);
%                 CpLiq1 = NitrousFluid.getLiquidCp(obj.temp,obj.vapourPressure);
%             end
%             gammaVap1 = NitrousFluid.getGasSpecificHeatRatio(obj.temp,obj.vapourPressure);
%             gammaLiq1 = NitrousFluid.getLiquidSpecificHeatRatio(obj.temp,obj.vapourPressure);
%             CvVap1 = CpVap1 / gammaVap1;
%             CvLiq1 = CpLiq1 / gammaLiq1;
%             %Commented stuff left for help making future improvement to
%             %unlinking temps of tank contents
%             
%             %Calculate dQ, the heat removed from liquid during this vapourization
% %             dQ = mVap * SaturatedNitrous.getLiquidEnthalpyOfVaporisation(obj.temp);
%             
%             %First let's remove the vapour mass that left the tank
%             dE = (-mDrained) * CpVap1 * obj.temp; %-Energy lost from vapour due to mass exiting tank (=mCpTo)
%             mVNew = obj.mVapour - mDrained; %New mass of vapour now some lost
%             
%             %Now let's calculate what happened to the liquid to vapourize
%             %mVap (in response to gas leaving)
%             %Calculate dT, the temp drop of the remaining liquid nitrous
%             %due to losing this heat (From vapourising). (dT = dQ/(m*C))
% %             dT = -dQ / ((obj.mLiquid) * CpLiq1);
% %             T2 = obj.temp +dT; %Calculate end temp of liquid
%             
%             mLNew = obj.mLiquid - mVap; %New mass of liquid
%             %Energy of contents of tank            
%             E1 = obj.mAir*CvAir*obj.temp + obj.mLiquid*CvLiq1*obj.temp + obj.mVapour*CvVap1*obj.temp + dE;
%             mVNew = mVNew + mVap; %Combine vapour masses (mVap and other vapour)
%             
%             %Uniform temperature of tank contents that gives same internal
%             %energy
%             T2 = E1 / (obj.mAir*CvAir + mLNew*CvLiq1 + mVNew*CvVap1);
%             
%             endState.temp = T2; %Update end temp of the end state
%             %Calculate mass vapourized in reaching end state, before liquid mass
%             %minus final liquid mass. Should equal mVap if mVap was correct
%             mVapourizedToReachEndState = (obj.mLiquid) - endState.mLiquid;
%             error = -(mVapourizedToReachEndState-mVap); %If positive then mVap was too big
%             %fprintf(['Vap: ',num2str(mVap),', err: ',num2str(error),' req Vap: ',num2str(mVapourizedToReachEndState),', T: ',num2str(T2),'\n']);
%         end
        
        %Numerically solves for physically correct volume of nitrous
        function updateVolNitrous(obj)            
            if obj.mAir == 0 %If no air then no need to solve anything complex
               obj.volNitrous = obj.tankTotalVolume;
                return;
            end
            if obj.mTotalNitrous == 0 %If no nitrous then no need to solve anything complex
               obj.volNitrous = 0;
               return;
            end

            initGuess= obj.volNitrous;
            if initGuess == 0
                initGuess = 0.1 * obj.tankTotalVolume;
            end
            solution = betterfzero(@pressureBoundaryError,initGuess,2*10^-7,obj.tankTotalVolume-2*10^-7,1*10^-2);
            obj.setVolNitrous(real(solution(1)));
            
            %Calculate the error in pressure difference at the boundary
            %between the air and the nitrous vapour for a given guess for
            %volNitrous
            function err = pressureBoundaryError(volNitrous)
                volNitrous = real(volNitrous);
                cloned = copy(obj);
                if volNitrous > obj.tankTotalVolume
                    volNitrous = obj.tankTotalVolume;
                end
                cloned.setVolNitrous(real(volNitrous));
                
                %Pressure at boundary between air and nitrous vapour, as
                %calulcated from P1 = P - rho*g*h
                P1 = cloned.vapourPressure - cloned.vapourDensity * GeometricNitrousAndAirTank.g * cloned.vapourHeight;
                %calculated from P2 = P + rho*g*h, assume air pressure
                %given (From PV=MRT) is the pressure at the mid-height of
                %the air, eg. the average pressure
                P2 = cloned.airPressure + cloned.airDensity * GeometricNitrousAndAirTank.g * (cloned.airHeight);
                err = P2-P1; %Should be zero if volNitrous given is physical
                %fprintf(['P vap: ',num2str(P1),' P air: ',num2str(P2),' Air h: ',num2str(cloned.airHeight),' airDens: ',num2str(cloned.airDensity),' volN: ',num2str(volNitrous),' \n']);
                %pause(0.1);
            end
        end
        
        %Define getter function for mass of liquid in tank
        function mLiq=getSaturatedLiquidMass(obj)
            vDensity = SaturatedNitrous.getVapourDensity(obj.temp);
            lDensity = SaturatedNitrous.getLiquidDensity(obj.temp);
            %From trivial re-arrangement of volume and mass summations
            mLiq = (obj.volNitrous - (obj.mTotalNitrous/vDensity)) / ((1/lDensity)-(1/vDensity));
        end
        
        %Define getter function for mass of vapour in tank
        function mVap=getSaturatedVapourMass(obj)
            vDensity = SaturatedNitrous.getVapourDensity(obj.temp);
            lDensity = SaturatedNitrous.getLiquidDensity(obj.temp);
            %From trivial re-arrangement of volume and mass summations
            mVap = (obj.volNitrous - (obj.mTotalNitrous/lDensity)) / ((1/vDensity)-(1/lDensity));
        end
    end
    
    methods
        %Create a nitrous tank object.
        %Args: (temp, mNitrousTotal, mAir,tankHeight,crossSectionalArea). Temp in KELVIN,
        %rest in standard SI units
        function obj = GeometricNitrousAndAirTank(temp,mNitrousTotal,mAir,tankHeight,crossSectionalArea,gasVentHoleHeight) %Constructor for object
            obj.tankTotalVolume = tankHeight*crossSectionalArea; %Set volume of tank
            obj.mAir = mAir;
            obj.mTotalNitrous = mNitrousTotal;
            obj.crossSectionArea = crossSectionalArea;
            obj.volNitrous = 0.5*obj.tankTotalVolume;
            obj.temp = temp; %Set temp of tank
            obj.gasVentHoleHeight = gasVentHoleHeight;
            obj.updateVolNitrous(); %Compute actual physical volume of nitrous
        end
        
        function E = getInternalEnergy(obj)
            E1 = obj.mAir * FluidType.AIR.getSpecificInternalEnergy(obj.temp,obj.airPressure);
            if obj.mAir == 0
               E1 = 0; %Stops equalling NaN
            end
            %PLiq = obj.getPressureAtHeight(0.5*obj.liquidHeight);
            E2 = obj.mLiquid * FluidType.NITROUS_LIQUID.getSpecificInternalEnergy(obj.temp,obj.vapourPressure);
            E3 = obj.mVapour * FluidType.NITROUS_GAS.getSpecificInternalEnergy(obj.temp,obj.vapourPressure);
            E = E1 + E2 + E3;
        end
        
        function E = getTankTotalEnthalpy(obj)
            E1 = obj.mAir * FluidType.AIR.getSpecificEnthalpy(obj.temp,obj.airPressure);
            if obj.mAir == 0
               E1 = 0; %Stops equalling NaN
            end
            %PLiq = obj.getPressureAtHeight(0.5*obj.liquidHeight);
            E2 = obj.mLiquid * FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.temp,obj.vapourPressure);
            E3 = obj.mVapour * FluidType.NITROUS_GAS.getSpecificEnthalpy(obj.temp,obj.vapourPressure);
            E = E1 + E2 + E3;
        end
        
        function drainAmountOfAir(obj,mDrain)
            %Number of small masses to calculate the drainage with
            numSegments = ceil(mDrain/GeometricNitrousAndAirTank.PREFERRED_SMALL_DRAIN_AMOUNT);
            dm = mDrain / numSegments; %Small amount of mass to drain with
            for i=1:(numSegments) %numSegments amount of times
                obj.drainSmallAmountOfAir(dm);
            end
        end
        
        function drainSmallAmountOfAir(obj,mDrain)
            h = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.temp,obj.airPressure);
            H = h*mDrain;
            obj.changeAirMassEnergy(-mDrain,-H);
        end
        
        %Sets the temperature of the tank to a given value ignoring the heat transfer required to do it, and then
        %updates the nitrous/air ratio of the tank (Which may alter the
        %temperature to be slightly different)
        function forceSetTemp(obj,temp)
            obj.temp = temp;
            obj.updateVolNitrous(); %May result in a change of the temperature
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
            obj.updateVolNitrous(); %Corrects volume of nitrous vs air (and changes temp if needed) in an energy conserving manner
            
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
        
        function Q = findHeatRateInputToKeepTempForDrainingLiquidRate(obj,mdotLiq)
            h = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.temp,obj.getPressureAtHeight(0));
            H = h*mdotLiq; %Energy being lost due to liquid
            dEdt = obj.findIntEnergyChangeRateForConstTemperatureWithFillRate(-mdotLiq); %dEdt to keep temp same
            Q = dEdt + H;
        end
        
        function Q = findHeatToMakeTemp(obj,tempDesired)
            clone = copy(obj);
            clone.temp = tempDesired;
            Q = clone.getInternalEnergy() - obj.getInternalEnergy();
        end
%         
%         function Q = findHeatToMakeTemp(obj,tempDesired)
%             c = NitrousFluid.getSaturatedHeatCapacity(obj.temp);
%             QEstimate = obj.mTotalNitrous * c * (tempDesired-obj.temp);
%             uncert = abs(0.3*QEstimate); %Assume +-30% of estimate
%             Q = real(betterfzero(@findErrQ,QEstimate,QEstimate-uncert,QEstimate+uncert));
%             
%             function errQ = findErrQ(Q)
%                 clone = copy(obj);
%                 clone.addHeat(real(Q));
%                 errQ = tempDesired - clone.temp;
%             end
%         end
        
        %Change the nitrous in the tank's mass and energy. Positive for
        %energy/mass input, and negative for energy/mass output. EIn is not
        %specific, but should take into account mass and type of each mass
        %in it's calculation
        function changeNitrousMassEnergy(obj,mIn,EIn)
            %Internal energy in the tank before
            clone = copy(obj);
            ECv1 = obj.getInternalEnergy();
            %disp("E Before: "+ECv1);
            ECv2 = ECv1 + EIn; %Internal energy of the tank after
%             disp("E change: "+EIn);
%             fprintf(['E1: ',num2str(ECv1),', E2: ',num2str(ECv2),'\n']);
%             fprintf(['M1: ',num2str(obj.mTotalNitrous),', M2: ',num2str(obj.mTotalNitrous + mIn),'\n']);
            
            obj.mTotalNitrous = obj.mTotalNitrous + mIn; %Nitrous mass of tank after
            clone.mTotalNitrous = obj.mTotalNitrous;
            %Solving for in saturated state max T until modelling can do
            %supercritical
            T2 = betterfzero(@errTemp,obj.temp,200,SaturatedNitrous.T_CRIT-2,2); %Solve for T that gives correct internal energy
            %Update tank temp and nitrous vol
            obj.temp = real(T2);
            obj.updateVolNitrous();
            
            function err = errTemp(T)       
               clone.temp = real(T);
                
               E2 = clone.getInternalEnergy();               
               err = E2 - ECv2;
%                fprintf(['T: ',num2str(T),', E: ',num2str(E2), ' Err: ',num2str(err),'\n']);
            end
        end
        
        %Change the air in the tank's mass and energy. Positive for
        %energy/mass input, and negative for energy/mass output. EIn is not
        %specific, but should take into account mass and type of each mass
        %in it's calculation
        function changeAirMassEnergy(obj,mIn,EIn)
            %Internal energy in the tank before
            before = copy(obj);
            ECv1 = before.getInternalEnergy();
            ECv2 = ECv1 + EIn; %Internal energy of the tank after
            obj.mAir = before.mAir + mIn; %Nitrous mass of tank after
            T2 = betterfzero(@errTemp,obj.temp,200,380,2); %Solve for T that gives correct internal energy
            %Update tank temp and nitrous vol
            obj.temp = real(T2);
            obj.updateVolNitrous();
            
            function err = errTemp(T)
               obj.temp = real(T);
               E2 = obj.getInternalEnergy();
               err = E2 - ECv2;
            end
        end
        
%         %Function to drain an amount of air out of the tank (Kg)
%         function drainSmallAmountOfAir(obj,mDrain)
%             before = copy(obj);
%             obj.mAir = obj.mAir - mDrain;
%             if obj.mAir < 0
%                 obj.mAir = 0;
%             end
%             
%             %Change temperature of tank based on energy conservation with
%             %small mass change of Cv
%             T = before.temp;
%             CvAir = GeometricNitrousAndAirTank.Cp_AIR - GeometricNitrousAndAirTank.R_AIR;
%             CvVap = NitrousFluid.getGasCv(T,before.vapourPressure); %Cv at vapour pressure is good enough
%             CvLiq = NitrousFluid.getLiquidCv(T,before.vapourPressure); %Cv at vapour pressure is good enough
%             %Energy in control volume (tank) before drain, assume all is
%             %thermal energy, assuming temperature is constant within Cv
%             %aswell
%             ECv1 = T * (before.mAir*CvAir + before.mVapour*CvVap + before.mLiquid*CvLiq);
%             dE = -mDrain * GeometricNitrousAndAirTank.Cp_AIR * T; %Energy lost due to lost air
%             ECv2 = ECv1 + dE; %New internal energy of tank
%             
%             %Solve for temp that gives final mVapour and mLiquid with
%             %correct internal energy
%             solution = betterfzero(@errFunc,obj.temp,1,378);
%             T2 = real(solution(1));
%             obj.temp = T2;
%             obj.updateVolNitrous(); %Corrects volume of nitrous vs air (and changes temp if needed) in an energy conserving manner
%             
%             function err = errFunc(endTemp)
%                 endTemp = real(endTemp);
%                 obj.temp = endTemp; %Changes masses of vapour and liquid within CV
%                 T2 = ECv2 / (obj.mAir*CvAir + obj.mVapour*CvVap + obj.mLiquid*CvLiq);
%                 err = T2 - endTemp;
%             end
%         end
        
        %Adds small amount of nitrous to the tank. EIn should be all energy input as a result of this mass addition, it should NOT be specific (NOT per Kg) 
        function addAmountOfNitrous(obj,mIn,EIn)
            %Number of small masses to calculate the drainage with
            numSegments = ceil(mIn/GeometricNitrousAndAirTank.PREFERRED_SMALL_DRAIN_AMOUNT);
            dm = mIn / numSegments; %Small amount of mass to drain with
            for i=1:(numSegments) %numSegments amount of times
                obj.addSmallAmountOfNitrous(dm,(EIn/numSegments));
            end
        end
        
        %Adds small amount of nitrous to the tank. EIn should be all energy input as a result of this mass addition, it should NOT be specific (NOT per Kg) 
        function addSmallAmountOfNitrous(obj,mIn,EIn)
            obj.changeNitrousMassEnergy(mIn,EIn);
        end
        
%         function addSmallAmountOfNitrous(obj,mIn,ho)            
%             before = copy(obj);
%             obj.mTotalNitrous = obj.mTotalNitrous + mIn;
%             %Change temperature of tank based on energy conservation with
%             %small mass change of Cv
%             T = before.temp;
%             CvAir = GeometricNitrousAndAirTank.Cp_AIR - GeometricNitrousAndAirTank.R_AIR;
%             CvVap = NitrousFluid.getGasCv(T,before.vapourPressure); %Cv at vapour pressure is good enough
%             CvLiq = NitrousFluid.getLiquidCv(T,before.vapourPressure); %Cv at vapour pressure is good enough
%             %Energy in control volume (tank) before drain, assume all is
%             %thermal energy, assuming temperature is constant within Cv
%             %aswell
%             ECv1 = T * (before.mAir*CvAir + before.mVapour*CvVap + before.mLiquid*CvLiq);
%             dE = mIn * ho;
%             ECv2 = ECv1 + dE;
%             
%             %Use guess for if was all vapour added as initial guess for
%             %temp
%             T2Approx = ECv2 / (obj.mAir*CvAir + (obj.mTotalNitrous-obj.mLiquid)*CvVap + obj.mLiquid*CvLiq);
%             
%             %Solve for temp that gives final mVapour and mLiquid with
%             %correct internal energy
%             solution = betterfzero(@errFunc,T2Approx,220,378);
%             T2 = real(solution(1));
%             obj.temp = T2;
%             obj.updateVolNitrous(); %Corrects volume of nitrous vs air (and changes temp if needed) in an energy conserving manner
%             
%             function err = errFunc(endTemp)
%                 endTemp = real(endTemp);
%                 obj.temp = endTemp; %Changes masses of vapour and liquid within CV
%                 %fprintf(['End temp guess: ',num2str(endTemp),' vapour pressure: ',num2str(obj.vapourPressure),'\n']);
%                 obj.updateVolNitrous();
%                 T2 = ECv2 / (obj.mAir*CvAir + obj.mVapour*CvVap + obj.mLiquid*CvLiq);
%                 err = T2 - endTemp;
%             end
%         end
%         
        %Re-calculates the state of the tank after draining an amount of
        %liquid from it.
        %Works by splitting amount into tiny amounts and then calling
        %drainSmallAmountOfLiquid(mDrained) until correct amount drained
        function drainAmountOfLiquid(obj,mDrained)
            %Number of small masses to calculate the drainage with
            numSegments = ceil(mDrained/GeometricNitrousAndAirTank.PREFERRED_SMALL_DRAIN_AMOUNT);
            dm = mDrained / numSegments; %Small amount of mass to drain with
            for i=1:(numSegments) %numSegments amount of times
                obj.drainSmallAmountOfLiquid(dm);
            end
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
            %Number of small masses to calculate the drainage with
            numSegments = ceil(mDrained/GeometricNitrousAndAirTank.PREFERRED_SMALL_DRAIN_AMOUNT);
            dm = mDrained / numSegments; %Small amount of mass to drain with
            for i=1:(numSegments) %numSegments amount of times
                obj.drainSmallAmountOfVapour(dm);
            end
        end
        
        function h = specificEnthalpyTopOfTank(obj)
            h = FluidType.NITROUS_GAS.getSpecificEnthalpy(obj.temp,obj.getPressureAtHeight(obj.gasVentHoleHeight));
        end
        
        function drainSmallAmountOfVapour(obj,mDrained)
           h = FluidType.NITROUS_GAS.getSpecificEnthalpy(obj.temp,obj.getPressureAtHeight(obj.gasVentHoleHeight));
           H = h*mDrained;
           obj.changeNitrousMassEnergy(-mDrained,-H);
        end
        
        function P = getPressureAtHeight(obj,h)
           P = obj.vapourPressure;
           nitrousHeight = obj.tankHeight - obj.airHeight;
           if h > nitrousHeight
              P = obj.airPressure; %If h lies within air, then uniform pressure is ok assumption
              return;
           end
           liqHeight = obj.liquidHeight; %Height of nitrous liquid phase
           if h < liqHeight
              P = P + obj.liquidDensity * GeometricNitrousAndAirTank.g * (liqHeight-h);
           else
              P = P - obj.vapourDensity * GeometricNitrousAndAirTank.g * (h-liqHeight);
           end
        end
        
%         % Re-calulates the state of the tank after draining a small amount
%         % of liquid from it. Has to be small as temp change approximated
%         % with current Hv, liquid mass and liquid specific heat capacity.
%         % Works by iteratively solving for mass that needs to have been vapourised
%         % For the system to stay saturated for a given reduction in mass. (Outflow
%         % of liquid)
%         % NOTE: Assumes boiling rate is infinitely high to keep tank in
%         % saturated state
%         function drainSmallAmountOfLiquid(obj, mDrained,mVapGuess)
%             %Use another nitrous tank object to help determine valid end
%             %state
%             endState = copy(obj);
%             endState.mTotalNitrous = obj.mTotalNitrous - mDrained; %Set end total mass
%             %Calculate mass vapourized by iteratively solving, to find
%             %correct end state
%             if ~exist('mVapGuess','var')
%                 mVapGuess = 1*10^-9;
%             end
%             mVap = mVapGuess; %Initial guess for amount of liquid vapourized
%             %Define tolerance for convergence of mVap, finer will be more
%             %accurate but also slower
%             maxError = max(mDrained / 10000, 1*10^-8);
%             %Solve via newton raphson method
%             solution = betterfzero(@errFunc,mVap,-mDrained*100-1*10^-6,mDrained*100+1*10^-6,maxError,2000);
%             mVap = real(solution(1));
%             errFunc(mVap);
%             %Uncomment to debug:
% %              mlStart = obj.mLiquid;
% %              T1 = obj.temp;
% %              mTotalBefore = obj.mTotal;
%             %Update the values of the actual tank
%             obj.mTotalNitrous = endState.mTotalNitrous;
%             obj.temp = endState.temp;
%             obj.updateVolNitrous(); %Corrects volume of nitrous vs air (and changes temp if needed) in an energy conserving manner
%             %Uncomment to debug:
%             %disp("Vap: "+mVap+", drain: "+mDrained+" mLiq: "+obj.mLiquid+", mLiqBefore: "+mlStart+" dMl: "+(obj.mLiquid-mlStart)+" T1: "+T1+" T2: "+obj.temp+" dMTotal: "+(obj.mTotal-mTotalBefore));
%             function err = errFunc(mVap)
%                 mVap = real(mVap);
%                err = calcVapourAmountGuess(obj,mDrained,endState,mVap); 
%             end
%         end
        
%         % Re-calulates the state of the tank after draining a small amount
%         % of vapour from it. Has to be small as temp change approximated
%         % with current Hv, liquid mass and liquid specific heat capacity.
%         % Works by iteratively solving for mass that needs to have been vapourised
%         % For the system to stay saturated for a given reduction in mass. (Outflow
%         % of vapour)
%         % NOTE: Assumes boiling rate is infinitely high to keep tank in
%         % saturated state
%         function drainSmallAmountOfVapour(obj, mDrained,mVapGuess)
%             %Use another nitrous tank object to help determine valid end
%             %state
%             endState = copy(obj);
%             endState.mTotalNitrous = obj.mTotalNitrous - mDrained; %Set end total mass
%             %Calculate mass vapourized by iteratively solving, to find
%             %correct end state
%             if ~exist('mVapGuess','var')
%                 mVapGuess = 1*10^-9;
%             end
%             mVap = mVapGuess; %Initial guess for amount of liquid vapourized
%             %Define tolerance for convergence of mVap, finer will be more
%             %accurate but also slower
%             maxError = max(mDrained / 100000, 1*10^-8);
%             %Solve via bracketed newton raphson method
%             solution = betterfzero(@errFunct,mVap,-mDrained*10^2-1*10^-6,mDrained*10^2+1*10^-6,maxError,1000);
%             mVap = real(solution(1));
%             errFunct(mVap); %Make sure value of end state correct
%             %Uncomment to debug:
% %              mlStart = obj.mLiquid;
% %              T1 = obj.temp;
% %              mTotalBefore = obj.mTotal;
%             %Update the values of the actual tank
%             obj.mTotalNitrous = endState.mTotalNitrous;
%             obj.temp = endState.temp;
%             obj.updateVolNitrous(); %Corrects volume of nitrous vs air (and changes temp if needed) in an energy conserving manner
%             %Uncomment to debug:
%             %disp("Vap: "+mVap+", drain: "+mDrained+" mLiq: "+obj.mLiquid+", mLiqBefore: "+mlStart+" dMl: "+(obj.mLiquid-mlStart)+" T1: "+T1+" T2: "+obj.temp+" dMTotal: "+(obj.mTotal-mTotalBefore));
%             function err = errFunct(mVap)
%                 mVap = real(mVap);
%                 err = calcVapourAmountGuessForDrainVapour(obj,mDrained,endState,mVap);
%             end
%         end
        
        function mTotalTank = get.mTotalTank(obj)
           mTotalTank = obj.mAir + obj.mTotalNitrous; 
        end
        
        %Define getter function for mass of liquid in tank
        function mLiq=get.mLiquid(obj)
            if ~obj.isSaturated() %If not at saturation conditions, then assume all nitrous is a gas
               mLiq = 0;
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
               mVap = obj.mTotalNitrous;
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
        
        %Define getter function for current pressure of the air at the top
        %of the tank
        function val=get.airPressure(obj)
            %P = rho*R*T
            val = obj.airDensity*GeometricNitrousAndAirTank.R_AIR*obj.temp;
        end
        
        %Define getter function for current air density
        function val=get.airDensity(obj)
            val = obj.mAir / obj.volAir;
        end
        
        %Define getter function for current air volume
        function val=get.volAir(obj)
            val = obj.tankTotalVolume - obj.volNitrous;
        end
        
        %Define getter function for current pressure exerted by vapour in
        %the tank
        function val=get.vapourPressure(obj)
            persistent lastVapourPressure;
            if ~obj.isSaturated()
                if obj.volNitrous == 0 %If no nitrous in the tank
                   val = obj.airPressure; %Just return air pressure
                   return;
                end
                
                if isempty(lastVapourPressure)
                    lastVapourPressure = obj.airPressure;
                end
               %Calculate pressure for given density with iteration...
               %solver = NonLinearEquationSolver({@errorNonSaturatedPressure},1*10^-7);
               %solution = solver.solve(obj.airPressure,[],0,SaturatedNitrous.getVapourPressure(obj.temp));
               solution = betterfzero(@errorNonSaturatedPressure,lastVapourPressure,1*10^-6,1*10^9,1*10^-7);
               val = real(solution(1)); %Solution is pressure in Pascal
               lastVapourPressure = val;
               return;
            end
            val = SaturatedNitrous.getVapourPressure(obj.temp);
            
            function err = errorNonSaturatedPressure(P)
                P = real(P);
                try
                    density = NitrousFluid.getGasDensity(obj.temp,P); %Density of the gas at the given temp and pressure
                catch err
                    disp("Failed P: "+P+" Failed temp: "+obj.temp);
                    rethrow(err);
                end
                if density < 0
                    density = 0;
                end
                err = density - obj.vapourDensity; %Error function is difference in density from what is required
                %fprintf(['Req density: ',num2str(obj.vapourDensity),', T: ',num2str(obj.temp),' P: ',num2str(P),' Density: ',num2str(density),' Err: ',num2str(err),'\n']);
                %fprintf(['P: ',num2str(P),'rho: ',num2str(density),' rho needed: ',num2str(obj.vapourDensity),'\n']);
                %pause(0.5);
            end
        end
        
        %Define getter function for current liquid density
        function val=get.liquidDensity(obj)
            val = SaturatedNitrous.getLiquidDensity(obj.temp);
        end
        
        %Define getter function for current vapour density
        function val=get.vapourDensity(obj)
            if ~obj.isSaturated() %If not saturated, then get density of nitrous gas
                %is just m/V
               val = obj.mTotalNitrous / obj.volNitrous;
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
            if ~obj.isSaturated() %If not saturated then vapour takes up entire tank
               vol = obj.volNitrous;
               return;
            end
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
        
        %Define getter function for current air height
        function val=get.airHeight(obj)
            val = obj.volAir / obj.crossSectionArea;
        end
        
        function val = get.isSaturated(obj)
            %If both >=0 then must be saturated, if one is negative then
            %not - instead assume nitrous is all vapour
            if true
                %TODO Fix without breaking vapour pressure calculations
                val = true;
                return;
            end
           val = obj.getSaturatedLiquidMass() >= 0 && obj.getSaturatedVapourMass() >= 0;
        end
        
        %Define getter function for the pressure at the base of the tank
        function val=get.pressureAtBaseOfTank(obj)
            %P1 = P + rho*g*h
            if obj.mTotalNitrous == 0
                val = obj.airPressure + obj.airDensity * GeometricNitrousAndAirTank.g * (obj.airHeight/2);
               return; 
            end
            val = obj.vapourPressure + obj.liquidDensity * GeometricNitrousAndAirTank.g * obj.liquidHeight;
        end
    end
end
