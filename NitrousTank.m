% Script to model the state of a nitrous tank, with given pertubations to mass etc.
% Depends on SaturatedNitrous.m
% By Eddie Brown

classdef NitrousTank < handle %Handle class so that passed by reference not value, makes behave sanely
    properties %Allowed to be modified
        temp; %Current temperature of the contents of the tank, IN KELVIN
    end
    properties (SetAccess=private) %Not allowed to be externally modified as may result in non-physical results
        mTotal; %Total mass of contents of the tank 
        tankTotalVol; %Tank's total volume
    end
    properties (Dependent) %Dependent properties
        mLiquid; %Current mass of liquid in the tank
        pressure; %The current pressure of the tank
        mVapour; %Current mass of vapour in the tank
        volVapour; %Current volume of the vapour in the tank
        volLiquid; %Current volume of the liquid in the tank
        liquidDensity; %Current density of the liquid in the tank
        vapourDensity; %Current density of the vapour in the tank
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
       PREFERRED_SMALL_DRAIN_AMOUNT = 0.0001; %0.1 grams
    end
    
    methods(Access=private)
        %Inner loop of the "drainSmallAmountOfLiquid"'s iterative solver
        %Uses guess for mass vapourized to calculate end state final
        %temperature and updates value for endState's temperature
        %accordingly.
        %This is then used to find the final liquid mass which is then used
        %to find the amount of liquid that must have vapourized
        %The error between the vapour amount guess and that calculated is
        %returned so that the iterative solver can find the solution
        %Is moved here as a separate function so can be called twice each
        %loop without code duplication (so gradient can be approximated for
        %newton raphson solver)
        function error = calcVapourAmountGuess(obj,mDrained,endState,mVap)
           %Calculate dQ, the heat removed from liquid during this vapourization
           dH = SaturatedNitrous.getLiquidEnthalpyOfVaporisation(obj.temp);
           dQ = mVap * dH;
           %Calculate dT, the temp drop of the remaining liquid nitrous
           %due to losing this heat. (dT = dQ/(m*C))
           dT = -dQ / ((obj.mLiquid - mDrained) * SaturatedNitrous.getLiquidIsobaricHeatCapacity(obj.temp));
           T2 = obj.temp +dT; %Calculate end temp
           endState.temp = T2; %Update end temp of the end state
           %Calculate mass vapourized in reaching end state, before liquid mass
           %minus final liquid mass. Should equal mVap if mVap was correct
           mVapourizedToReachEndState = (obj.mLiquid - mDrained) - endState.mLiquid;
           error = -(mVapourizedToReachEndState-mVap); %If positive then mVap was too big 
        end
        
        %Inner loop of the "drainSmallAmountOfValpour"'s iterative solver
        %Uses guess for mass vapourized to calculate end state final
        %temperature and updates value for endState's temperature
        %accordingly.
        %This is then used to find the final liquid mass which is then used
        %to find the amount of liquid that must have vapourized
        %The error between the vapour amount guess and that calculated is
        %returned so that the iterative solver can find the solution
        %Is moved here as a separate function so can be called twice each
        %loop without code duplication (so gradient can be approximated for
        %newton raphson solver)
        function error = calcVapourAmountGuessForDrainVapour(obj,endState,mVap)
           %Calculate dQ, the heat removed from liquid during this vapourization
           dQ = mVap * SaturatedNitrous.getLiquidEnthalpyOfVaporisation(obj.temp);
           %Calculate dT, the temp drop of the remaining liquid nitrous
           %due to losing this heat. (dT = dQ/(m*C))
           dT = -dQ / ((obj.mLiquid) * SaturatedNitrous.getLiquidIsobaricHeatCapacity(obj.temp));
           T2 = obj.temp +dT; %Calculate end temp
           endState.temp = T2; %Update end temp of the end state
           %Calculate mass vapourized in reaching end state, before liquid mass
           %minus final liquid mass. Should equal mVap if mVap was correct
           mVapourizedToReachEndState = (obj.mLiquid) - endState.mLiquid;
           error = -(mVapourizedToReachEndState-mVap); %If positive then mVap was too big 
        end
    end
    
    methods
        %Create a nitrous tank object.
        %Args: (temp, mLiquid, tankTotalVol). Temp in KELVIN,
        %rest in SI units
        function obj = NitrousTank(temp,mLiquid,tankTotalVol) %Constructor for object
            obj.temp = temp; %Set temp of tank
            obj.tankTotalVol = tankTotalVol; %Set volume of tank
            volLiqInit = mLiquid / obj.liquidDensity; %Calculate vol of liquid
            volVapInit = tankTotalVol - volLiqInit; %Calculate vol of vapour
            obj.mTotal = mLiquid + volVapInit*obj.vapourDensity; %Set total mass of tank from calculated variables
        end
        
        %Re-calculates the state of the tank after draining an amount of
        %liquid from it.
        %Works by splitting amount into tiny amounts and then calling
        %drainSmallAmountOfLiquid(mDrained) until correct amount drained
        function drainAmountOfLiquid(obj,mDrained)
            %Number of small masses to calculate the drainage with
            numSegments = ceil(mDrained/NitrousTank.PREFERRED_SMALL_DRAIN_AMOUNT);
            dm = mDrained / numSegments; %Small amount of mass to drain with
            totalDrained = 0;
            mLiqStart = obj.mLiquid;
            for i=1:(numSegments) %numSegments amount of times
                obj.drainSmallAmountOfLiquid(dm);
                totalDrained = totalDrained + dm;
            end
            if(obj.mLiquid > mLiqStart-mDrained)
               error("Mass drainage solver failed! Drained mass was not all liquid!"); 
            end
        end
        
        %Re-calculates the state of the tank after draining an amount of
        %vapour from it.
        %Works by splitting amount into tiny amounts and then calling
        %drainSmallAmountOfVapour(mDrained) until correct amount drained
        function drainAmountOfVapour(obj,mDrained)
            %Number of small masses to calculate the drainage with
            numSegments = ceil(mDrained/NitrousTank.PREFERRED_SMALL_DRAIN_AMOUNT);
            dm = mDrained / numSegments; %Small amount of mass to drain with
            %totalDrained = 0;
            for i=1:(numSegments) %numSegments amount of times
                obj.drainSmallAmountOfVapour(dm);
               % totalDrained = totalDrained + dm;
            end
        end
        
        % Re-calulates the state of the tank after draining a small amount
        % of liquid from it. Has to be small as temp change approximated
        % with current Hv, liquid mass and liquid specific heat capacity.
        % Works by iteratively solving for mass that needs to have been vapourised
        % For the system to stay saturated for a given reduction in mass. (Outflow
        % of liquid)
        % NOTE: Assumes boiling rate is infinitely high to keep tank in
        % saturated state
        function drainSmallAmountOfLiquid(obj, mDrained)
            %Use another nitrous tank object to help determine valid end
            %state
            endState = NitrousTank(obj.temp,obj.mLiquid,obj.tankTotalVol);
            endState.mTotal = obj.mTotal - mDrained; %Set end total mass
            %Calculate mass vapourized by iteratively solving, to find
            %correct end state
            error = 999999;
            mVap = 1*10^-6; %Initial guess for amount of liquid vapourized
            %Define tolerance for convergence of mVap, finer will be more
            %accurate but also slower
            maxError = min(mDrained / 10000, 1*10^-8);
            iter = 0;
            %Solve via newton raphson method
            while (abs(error) > maxError && iter < 10000) %Error of 0.01 milligram allowed, max 10k iterations
                iter = iter + 1;
                errorSmallPertub = calcVapourAmountGuess(obj,mDrained,endState,mVap+1*10^-7); %Value of error function with small pertubation from guess
                error = calcVapourAmountGuess(obj,mDrained,endState,mVap); %See function for what it does, value of error function, updates endState's final temp within function
                gradientApprox = (errorSmallPertub-error) / (1*10^-7); %dy/dx
                mVap = mVap-(error/gradientApprox); %Adjust mVap using newton raphson method
            end
            if (abs(error) > maxError)
                warning("Mass drainage iterative solver failed to converge! Solution will be invalid! (Maybe tank is empty?). Liquid pre-small-drain: "+obj.mLiquid);
            end
            if (mVap < 0) %Gas condensing not physical
                error("NEGATIVE mVap");
            end
            if (endState.mLiquid > obj.mLiquid - mDrained)
                error("Mass drainage solver failed! Drained mass was not liquid!");
            end
            %Uncomment to debug:
%              mlStart = obj.mLiquid;
%              T1 = obj.temp;
%              mTotalBefore = obj.mTotal;
            %Update the values of the actual tank
            obj.mTotal = endState.mTotal;
            obj.temp = endState.temp;
            %Uncomment to debug:
            %disp("Vap: "+mVap+", drain: "+mDrained+" mLiq: "+obj.mLiquid+", mLiqBefore: "+mlStart+" dMl: "+(obj.mLiquid-mlStart)+" T1: "+T1+" T2: "+obj.temp+" dMTotal: "+(obj.mTotal-mTotalBefore));
        end
        
        % Re-calulates the state of the tank after draining a small amount
        % of vapour from it. Has to be small as temp change approximated
        % with current Hv, liquid mass and liquid specific heat capacity.
        % Works by iteratively solving for mass that needs to have been vapourised
        % For the system to stay saturated for a given reduction in mass. (Outflow
        % of vapour)
        % NOTE: Assumes boiling rate is infinitely high to keep tank in
        % saturated state
        function drainSmallAmountOfVapour(obj, mDrained)
            %Use another nitrous tank object to help determine valid end
            %state
            endState = NitrousTank(obj.temp,obj.mLiquid,obj.tankTotalVol);
            endState.mTotal = obj.mTotal - mDrained; %Set end total mass
            %Calculate mass vapourized by iteratively solving, to find
            %correct end state
            error = 999999;
            mVap = 1*10^-6; %Initial guess for amount of liquid vapourized
            %Define tolerance for convergence of mVap, finer will be more
            %accurate but also slower
            maxError = min(mDrained / 10000, 1*10^-8);
            iter = 0;
            %Solve via newton raphson method
            while (abs(error) > maxError && iter < 10000) %Error of 0.01 milligram allowed, max 10k iterations
                iter = iter + 1;
                errorSmallPertub = calcVapourAmountGuessForDrainVapour(obj,endState,mVap+1*10^-7); %Value of error function with small pertubation from guess
                error = calcVapourAmountGuessForDrainVapour(obj,endState,mVap); %See function for what it does, value of error function, updates endState's final temp within function
                gradientApprox = (errorSmallPertub-error) / (1*10^-7); %dy/dx
                mVap = mVap-(error/gradientApprox); %Adjust mVap using newton raphson method
            end
            if (abs(error) > maxError)
                warning("Mass drainage iterative solver failed to converge! Solution will be invalid! (Maybe tank is empty?). Liquid pre-small-drain: "+obj.mLiquid);
            end
            if (mVap < 0) %Gas condensing not physical
                error("NEGATIVE mVap");
            end
            if (endState.mLiquid > obj.mLiquid)
                error("Mass drainage solver failed! Increase in liquid within tank!");
            end
            %Uncomment to debug:
%              mlStart = obj.mLiquid;
%              T1 = obj.temp;
%              mTotalBefore = obj.mTotal;
            %Update the values of the actual tank
            obj.mTotal = endState.mTotal;
            obj.temp = endState.temp;
            %Uncomment to debug:
            %disp("Vap: "+mVap+", drain: "+mDrained+" mLiq: "+obj.mLiquid+", mLiqBefore: "+mlStart+" dMl: "+(obj.mLiquid-mlStart)+" T1: "+T1+" T2: "+obj.temp+" dMTotal: "+(obj.mTotal-mTotalBefore));
        end
        
        %Define getter function for mass of liquid in tank
        function mLiq=get.mLiquid(obj)
            %From trivial re-arrangement of volume and mass summations
            mLiq = (obj.tankTotalVol - (obj.mTotal/obj.vapourDensity)) / ((1/obj.liquidDensity)-(1/obj.vapourDensity));
            if (mLiq < 0)
               warning("Liquid mass negative, no physical saturated solution can exist for these conditions"); 
            end
        end
        
        %Define getter function for mass of vapour in tank
        function mVap=get.mVapour(obj)
            %From trivial re-arrangement of volume and mass summations
            mVap = (obj.tankTotalVol - (obj.mTotal/obj.liquidDensity)) / ((1/obj.vapourDensity)-(1/obj.liquidDensity));
            if (mVap < 0)
               warning("Vapour mass negative, no physical saturated solution can exist for these conditions"); 
            end
        end
        
        %Define getter function for current tank pressure
        function val=get.pressure(obj)
            val = SaturatedNitrous.getVapourPressure(obj.temp);
        end
        
        %Define getter function for current liquid density
        function val=get.liquidDensity(obj)
            val = SaturatedNitrous.getLiquidDensity(obj.temp);
        end
        
        %Define getter function for current vapour density
        function val=get.vapourDensity(obj)
            val = SaturatedNitrous.getVapourDensity(obj.temp);
        end
        
        %Define getter function for current liquid volume
        function vol=get.volLiquid(obj)
           vol = obj.mLiquid / obj.liquidDensity; %V = m/rho
        end
        
        %Define getter function for current vapour volume
        function vol=get.volVapour(obj)
           vol = obj.mVapour / obj.vapourDensity; %V = m/rho
        end
    end
end
