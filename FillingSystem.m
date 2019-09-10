classdef FillingSystem  < matlab.mixin.Copyable
   properties
      internalTank;
      externalTank;
      pipeBetweenTanks;
      ventPipe;
      fillValveOpenAmt = 1;
      ventValveOpenAmt = 0;
      ambientTemp = 25+FillingSystem.C_TO_K;
      QExt = 0; %Power (J/s) of heat being applied to external tank system by heating/cooling
      internalTankSurfaceArea;
      heatTransferCoeffInternalTank;    
      externalTankSystem;
   end
   properties (Constant)
      C_TO_K = 273.15; 
      ATM_PRESSURE = 101325;
   end
   methods(Access = protected)
       % Override copyElement method:
       function cpObj = copyElement(obj)
           % Make a shallow copy of all primitive properties
           cpObj = copyElement@matlab.mixin.Copyable(obj);
           % Make a deep copy of the necessary objects
           cpObj.internalTank = copy(obj.internalTank);
           cpObj.externalTank = copy(obj.externalTank);
           cpObj.pipeBetweenTanks = copy(obj.pipeBetweenTanks);
           cpObj.ventPipe = copy(obj.ventPipe);
           cpObj.externalTankSystem = copy(obj.externalTankSystem);
       end
   end
   methods
       function obj = FillingSystem(ambientTemp,internalTankHeight,internalTankCrossSectionA...
               ,internalVentHoleHeight,initialInternalNitrousMass,initialExternalNitrousMass,...
               externalTankHeight,externalTankCrossSectionA,internalTankSurfaceArea,...
               heatTransferCoeffInternalTank,externalTankSurfaceArea,...
               heatTransferCoeffExternalTankWithFluid,SHCFluid,externalTankFluidSurfaceArea,...
               externalTankFluidHeatTransferCoeffWithSurroundings,mFluid,initialInternalTankTemp,...
               initialExternalTankTemp)
            obj.ambientTemp = ambientTemp;
            obj.internalTank = GeometricNitrousTank(initialInternalTankTemp,initialInternalNitrousMass,internalTankHeight,internalTankCrossSectionA,internalVentHoleHeight); 
            obj.externalTank = GeometricNitrousTank(initialExternalTankTemp,initialExternalNitrousMass,externalTankHeight,externalTankCrossSectionA,externalTankHeight);
            obj.ventPipe = SaturatedPipeValvePipeFast('preBakedData/saturatedGasPipeValveFlowRates.mat');
            obj.pipeBetweenTanks = SaturatedPipeValvePipeFast('preBakedData/saturatedLiquidPipeValveFlowRates.mat');
            obj.internalTankSurfaceArea = internalTankSurfaceArea;
            obj.heatTransferCoeffInternalTank = heatTransferCoeffInternalTank;
            obj.externalTankSystem = ExternalTankSystem(obj.externalTank,externalTankSurfaceArea,...
                heatTransferCoeffExternalTankWithFluid,ambientTemp, SHCFluid, ...
               externalTankFluidSurfaceArea,externalTankFluidHeatTransferCoeffWithSurroundings,...
               mFluid);
       end
       
       function [fillValveOpenAmt,ventValveOpenAmt,TFluidExtReq,Q,mdotBetweenTanks,mdotVent,QInclFromEnv] = findControlPointForConditions(obj,mdotFillRate,internalTempChangeRate,externalTempChangeRate,externalFluidTempChangeRate)
           [mdotBetweenTanks,mdotVent,TFluidExtReq,Q,QInclFromEnv] = obj.findPointForConditions(mdotFillRate,internalTempChangeRate,externalTempChangeRate,externalFluidTempChangeRate);
           try
               fillValveOpenAmt = betterfzero(@getMDotBetweenTanksErr,0.5,0,1,1e-7,1000,0.5*10^-16,true);
           catch excep
               disp("Error finding fill valve open amt for mass flow "+mdotBetweenTanks);
               rethrow excep;
           end
           try
               ventValveOpenAmt = betterfzero(@getMDotVentErr,0.5,0,1,1e-7,1000,0.5*10^-16,true);
           catch excep
               disp("Error finding vent valve open amt for mass flow "+mdotVent);
               rethrow excep;
           end
           
           function mdotBetweenTanksErr = getMDotBetweenTanksErr(fillValveOpenAmt)
               mdotBetweenTanksErr = mdotBetweenTanks - obj.getFlowBetweenTanksIfFillValveOpenAmtWasThis(fillValveOpenAmt);
           end
           
           function mdotVentErr = getMDotVentErr(ventValveOpenAmt)
               mdotVentErr = mdotVent - obj.getFlowOutIfVentValveOpenAmtWasThis(ventValveOpenAmt);
           end
       end
       
       function [fillValveOpenAmt,ventValveOpenAmt,TFluidExtReq,Q,mdotBetweenTanks,mdotVent,QInclFromEnv] = findControlEquilibriumPointForFillRate(obj,mdotFillRate)
           [mdotBetweenTanks,mdotVent,TFluidExtReq,Q,QInclFromEnv] = obj.findEquilibriumPointForFillRate(mdotFillRate);
           %Find fillValveOpenAmt that gives mdotBetweenTanks
           
%            disp("Between tanks: "+mdotBetweenTanks);
%            disp("Vent: "+mdotVent);
%            disp("T Ext fluid: "+TFluidExtReq);
%            disp("Q: "+Q);
           
            try
           fillValveOpenAmt = betterfzero(@getMDotBetweenTanksErr,0.5,0,1,1e-7,1000,0.5*10^-16,true);
            catch excep
                disp("Error finding fill valve open amt for mass flow "+mdotBetweenTanks);
                rethrow excep;
            end
            try
                ventValveOpenAmt = betterfzero(@getMDotVentErr,0.5,0,1,1e-7,1000,0.5*10^-16,true);
            catch excep
                disp("Error finding vent valve open amt for mass flow "+mdotVent);
                rethrow excep;
            end
           
           function mdotBetweenTanksErr = getMDotBetweenTanksErr(fillValveOpenAmt)
               mdotBetweenTanksErr = mdotBetweenTanks - obj.getFlowBetweenTanksIfFillValveOpenAmtWasThis(fillValveOpenAmt);
           end
           
           function mdotVentErr = getMDotVentErr(ventValveOpenAmt)
               mdotVentErr = mdotVent - obj.getFlowOutIfVentValveOpenAmtWasThis(ventValveOpenAmt);
           end
       end
       
       function [mdotBetweenTanks,mdotVent,TFluidExtReq,Q,QIncludingFromEnv] = findPointForConditions(obj,mdotFillRate,internalTempChangeRate,externalTempChangeRate,externalFluidTempChangeRate)
           mdotVent = obj.calcMDotOutThisInternalTankReqInclTempChangeRate(mdotFillRate,internalTempChangeRate);
           mdotBetweenTanks = mdotVent + mdotFillRate;
           %Calculate Q
           QTankReq = obj.externalTankSystem.externalTank.findHeatRateInputForTempChangeRateWithLiquidDrainingRate(externalTempChangeRate,mdotFillRate);
           TFluidExtReq = obj.externalTankSystem.calcFluidTempRequiredForHeatFluxToTank(QTankReq);
           [Q,QIncludingFromEnv] = obj.externalTankSystem.calcExternalHeatFluxReqForFluidTempChangeRate(TFluidExtReq,obj.ambientTemp,externalFluidTempChangeRate);
       end
       
       function [mdotBetweenTanks,mdotVent,TFluidExtReq,Q,QIncludingFromEnv] = findEquilibriumPointForFillRate(obj,mdotFillRate)
           mdotVent = obj.calcMDotOutThisInternalTankReq(mdotFillRate);
           mdotBetweenTanks = mdotVent + mdotFillRate;
           %Calculate Q
           QTankReq = obj.externalTankSystem.externalTank.findHeatRateInputToKeepTempForDrainingLiquidRate(mdotFillRate);
           TFluidExtReq = obj.externalTankSystem.calcFluidTempRequiredForHeatFluxToTank(QTankReq);
           [Q,QIncludingFromEnv] = obj.externalTankSystem.calcExternalHeatFluxReqToKeepFluidAtTemp(TFluidExtReq,obj.ambientTemp);
       end
       
       function mdotOut = calcMDotOutThisInternalTankReq(obj,mdotFillRate)
            heatTrans = obj.calcCurrentHeatTransferRateAmbientToInternalTank();
            mdotOut = obj.calcMDotOutInternalTankReq(heatTrans,mdotFillRate,obj.internalTank,obj.externalTankSystem.externalTank.temp);
       end
       
       function mdotOut = calcMDotOutThisInternalTankReqInclTempChangeRate(obj,mdotFillRate,internalTankTemperatureChangeRate)
            heatTrans = obj.calcCurrentHeatTransferRateAmbientToInternalTank();
            mdotOut = obj.calcMDotOutInternalTankReqInclTempChangeRate(heatTrans,mdotFillRate,obj.internalTank,obj.externalTankSystem.externalTank.temp,internalTankTemperatureChangeRate);
       end
       
       function mdotOut = calcMDotOutInternalTankReq(obj,HeatTransferIntoTank,mdotFillRate,tank,extTankT)
           hFromExt = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,extTankT,FluidProperty.VAPOR_QUALITY,0);
           hInt = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,tank.temp,FluidProperty.VAPOR_QUALITY,1);
           dEdt = tank.findIntEnergyChangeRateForConstTemperatureWithFillRate(mdotFillRate);
           mdotOut = (HeatTransferIntoTank - dEdt + mdotFillRate.*hFromExt)/(hInt - hFromExt);
       end
       
       function mdotOut = calcMDotOutInternalTankReqInclTempChangeRate(obj,HeatTransferIntoTank,mdotFillRate,tank,extTankT,internalTankTemperatureChangeRate)
           hFromExt = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,extTankT,FluidProperty.VAPOR_QUALITY,0);
           hInt = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,tank.temp,FluidProperty.VAPOR_QUALITY,1);
           dEdt = tank.findIntEnergyChangeRateForTempChangeRateWithFillRate(internalTankTemperatureChangeRate,mdotFillRate);
           mdotOut = (HeatTransferIntoTank - dEdt + mdotFillRate.*hFromExt)/(hInt - hFromExt);
       end
       
       function xDot = getSystemChangeRateVector(obj)
            x1 = obj.getStateVector();
            obj2 = copy(obj);
            dt = 0.01;
            obj2.advanceSystemBySmallTimeIncrem(dt);
            x2 = obj2.getStateVector();
            xDot = (x2-x1)./dt;
       end
       
       function QRate = calcCurrentHeatTransferRateAmbientToInternalTank(obj)
            QRate = obj.internalTankSurfaceArea .* obj.heatTransferCoeffInternalTank .* (obj.ambientTemp - obj.internalTank.temp);
       end
       
       function mdotBetween = getFlowBetweenTanksIfFillValveOpenAmtWasThis(obj,fillvalveOpenAmt)
            mdotBetween = obj.pipeBetweenTanks.getMassFlow(obj.externalTank.vapourPressure,obj.internalTank.vapourPressure,fillvalveOpenAmt);
       end
       
       function mdotVent = getFlowOutIfVentValveOpenAmtWasThis(obj,ventValveOpenAmt)
           mdotVent = obj.ventPipe.getMassFlow(obj.internalTank.vapourPressure,FillingSystem.ATM_PRESSURE,ventValveOpenAmt);
       end
       
       function advanceSystemBySmallTimeIncrem(obj,dt)
            mdotFill = obj.pipeBetweenTanks.getMassFlow(obj.externalTank.vapourPressure,obj.internalTank.vapourPressure,obj.fillValveOpenAmt);
            mdotVent = obj.ventPipe.getMassFlow(obj.internalTank.vapourPressure,FillingSystem.ATM_PRESSURE,obj.ventValveOpenAmt);
            EFlowBetweenTanks = dt.* mdotFill .* FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.externalTank.temp,obj.exernalTank.vapourPressure);
            obj.externalTank.changeNitrousMassEnergy(-mdotFill.*dt,-EFlowBetweenTanks);
            obj.internalTank.changeNitrousMassEnergy(mdotFill.*dt,EFlowBetweenTanks);
            if(obj.internalTank.liquidHeight >= obj.internalTank.gasVentHoleHeight)
               hIntVenting = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.internalTank.temp,obj.internalTank.vapourPressure);
            else
               hIntVenting = FluidType.NITROUS_GAS.getSpecificEnthalpy(obj.internalTank.temp,obj.internalTank.vapourPressure); 
            end
            EFlowOutOfRunTank = dt.* mdotVent .* hIntVenting;
            obj.internalTank.changeNitrousMassEnergy(-mdotVent.*dt,-EFlowOutOfRunTank);
            %HEAT TRANSFERS
            QInternal = dt .* obj.calcCurrentHeatTransferRateAmbientToInternalTank();
            obj.internalTank.addHeat(QInternal);
            obj.externalTankSystem.doHeatTransferForTimeStep(dt,obj.ambientTemp,obj.QExt);
       end
       
       function obj2 = getCopyInNewState(obj,x,u)
           obj2 = copy(obj);
           obj2.ambientTemp = x(1);
           obj2.externalTank.forceSetTemp(x(2));
           obj2.externalTank.forceSetNitrousMass(x(3));
           obj2.internalTank.forceSetTemp(x(4));
           obj2.internalTank.forceSetNitrousMass(x(5));
           obj2.externalTankSystem.TFluid = x(6);
           obj2.fillValveOpenAmt = u(1);
           obj2.ventValveOpenAmt = u(2);
           obj2.QExt = u(3);
       end
       
       function x = getStateVector(obj)
           x = [obj.ambientTemp;
               obj.externalTank.temp;
               obj.externalTank.mTotalNitrous;
               obj.internalTank.temp;
               obj.internalTank.mTotalNitrous;
               obj.externalTankSystem.TFluid];
       end
       
       function u = getControlVector(obj)
           u = [obj.fillValveOpenAmt;
               obj.ventValveOpenAmt;
               obj.QExt];
       end
   end
end