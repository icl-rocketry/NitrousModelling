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
           cpObj.externalTankSystem.externalTank = cpObj.externalTank;
       end
   end
   methods
       function obj = FillingSystem(ambientTemp,internalTankHeight,internalTankCrossSectionA...
               ,internalVentHoleHeight,initialInternalNitrousMass,initialExternalNitrousMass,...
               externalTankHeight,externalTankCrossSectionA,internalTankSurfaceArea,...
               heatTransferCoeffInternalTank,externalTankSurfaceArea,...
               heatTransferCoeffExternalTankWithFluid,SHCFluid,externalTankFluidSurfaceArea,...
               externalTankFluidHeatTransferCoeffWithSurroundings,mFluid,initialInternalTankTemp,...
               initialExternalTankTemp,initialExtFluidTemp)
            obj.ambientTemp = ambientTemp;
            obj.internalTank = GeometricNitrousTank(initialInternalTankTemp,initialInternalNitrousMass,internalTankHeight,internalTankCrossSectionA,internalVentHoleHeight); 
            obj.externalTank = GeometricNitrousTank(initialExternalTankTemp,initialExternalNitrousMass,externalTankHeight,externalTankCrossSectionA,externalTankHeight);
            obj.ventPipe = SaturatedPipeValvePipeFast('preBakedData/saturatedGasPipeValveFlowRates.mat','preBakedData/saturatedGasPipeValveFlowRatesInverted.mat');
            obj.pipeBetweenTanks = SaturatedPipeValvePipeFast('preBakedData/saturatedLiquidPipeValveFlowRates.mat','preBakedData/saturatedLiquidPipeValveFlowRatesInverted.mat');
            obj.internalTankSurfaceArea = internalTankSurfaceArea;
            obj.heatTransferCoeffInternalTank = heatTransferCoeffInternalTank;
            obj.externalTankSystem = ExternalTankSystem(obj.externalTank,externalTankSurfaceArea,...
                heatTransferCoeffExternalTankWithFluid,initialExtFluidTemp, SHCFluid, ...
               externalTankFluidSurfaceArea,externalTankFluidHeatTransferCoeffWithSurroundings,...
               mFluid);
       end
       
       function fillValveOpenAmt = findFillValvePositionForFlow(obj,mdotBetweenTanks)
           fillValveOpenAmt = 0.5;
           if(obj.pipeBetweenTanks.hasInvertedData)
                fillValveOpenAmt = obj.pipeBetweenTanks.getValvePos(obj.externalTank.vapourPressure,obj.internalTank.vapourPressure,mdotBetweenTanks);
               return;
           end
           try
               fillValveOpenAmt = betterfzero(@getMDotBetweenTanksErr,fillValveOpenAmt,0,1,1e-7,1000,0.5*10^-16,true);
           catch excep
               disp("Error finding fill valve open amt for mass flow "+mdotBetweenTanks);
               rethrow excep;
           end
           function mdotBetweenTanksErr = getMDotBetweenTanksErr(fillValveOpenAmt)
               mdotBetweenTanksErr = mdotBetweenTanks - obj.getFlowBetweenTanksIfFillValveOpenAmtWasThis(fillValveOpenAmt);
           end
       end
       
       function ventValveOpenAmt = findVentValvePositionForFlow(obj,mdotVent)
           ventValveOpenAmt = 0.5;
           if(obj.ventPipe.hasInvertedData)
                ventValveOpenAmt = obj.ventPipe.getValvePos(obj.internalTank.vapourPressure,FillingSystem.ATM_PRESSURE,mdotVent);
               return;
           end
           try
               ventValveOpenAmt = betterfzero(@getMDotVentErr,ventValveOpenAmt,0,1,1e-7,1000,0.5*10^-16,true);
           catch except
               disp("Error finding vent valve open amt for mass flow "+mdotVent);
               disp(except);
               rethrow except;
           end
           function mdotVentErr = getMDotVentErr(ventValveOpenAmt)
               mdotVentErr = mdotVent - obj.getFlowOutIfVentValveOpenAmtWasThis(ventValveOpenAmt);
           end
       end
       
       function [fillValveOpenAmt,dTintdt,TFluidExtReq,Q,QInclFromEnv] = calcFillingPointNoVentExternalTankEquilibrium(obj,mdotBetweenTanks)
           fillValveOpenAmt = obj.findFillValvePositionForFlow(mdotBetweenTanks);
           %Calculate Q
           QTankReq = obj.externalTankSystem.externalTank.findHeatRateInputForTempChangeRateWithLiquidDrainingRate(0,mdotBetweenTanks);
           TFluidExtReq = obj.externalTankSystem.calcFluidTempRequiredForHeatFluxToTank(QTankReq);
           [Q,QInclFromEnv] = obj.externalTankSystem.calcExternalHeatFluxReqForFluidTempChangeRate(TFluidExtReq,obj.ambientTemp,0);
           systemCopy = copy(obj);
           systemCopy.fillValveOpenAmt = fillValveOpenAmt;
           systemCopy.ventValveOpenAmt = 0;
           systemCopy.QExt = Q;
           systemCopy.externalTankSystem.TFluid = TFluidExtReq;
           xdot = systemCopy.getSystemChangeRateVector();
           dTintdt = xdot(4);
       end
       
       function [fillValveOpenAmt,dTintdt,dTextdt,dTextFdt] = calcEffectOfFillingNoVentWithHeat(obj,mdotBetweenTanks,QExt)
           fillValveOpenAmt = obj.findFillValvePositionForFlow(mdotBetweenTanks);
           systemCopy = copy(obj);
           systemCopy.fillValveOpenAmt = fillValveOpenAmt;
           systemCopy.ventValveOpenAmt = 0;
           systemCopy.QExt = QExt;
           xdot = systemCopy.getSystemChangeRateVector();
           dTintdt = xdot(4);
           dTextdt = xdot(2);
           dTextFdt = xdot(6);
       end
       
       function xdot = calcEffectOfControlOnSystem(obj,fillValveOpenAmt,ventValveOpenAmt,Q)
           systemCopy = copy(obj);
           systemCopy.fillValveOpenAmt = fillValveOpenAmt;
           systemCopy.ventValveOpenAmt = ventValveOpenAmt;
           systemCopy.QExt = Q;
           xdot = systemCopy.getSystemChangeRateVector();
       end
       
       function xdot = calcEffectOfControlOnSystemWithDiffExtFluidT(obj,fillValveOpenAmt,ventValveOpenAmt,Q,TExtFluid)
           systemCopy = copy(obj);
           systemCopy.externalTankSystem.TFluid = TExtFluid;
           systemCopy.fillValveOpenAmt = fillValveOpenAmt;
           systemCopy.ventValveOpenAmt = ventValveOpenAmt;
           systemCopy.QExt = Q;
           xdot = systemCopy.getSystemChangeRateVector();
           disp("TFluid: "+TExtFluid);
           disp(xdot);
       end
       
       function [dTintdt] = calcEffectOfFullyVentingNoFill(obj)
           systemCopy = copy(obj);
           systemCopy.fillValveOpenAmt = 0;
           systemCopy.ventValveOpenAmt = 1;
           xdot = systemCopy.getSystemChangeRateVector();
           dTintdt = xdot(4);
       end
       
       function [fillValveOpenAmt,dTintdt,dTextdt,dTextFdt] = calcEffectOfFillingNoVentOrHeat(obj,mdotBetweenTanks)
           [fillValveOpenAmt,dTintdt,dTextdt,dTextFdt] = obj.calcEffectOfFillingNoVentWithHeat(mdotBetweenTanks,0);
       end
       
       function [fillValveOpenAmt,ventValveOpenAmt,TFluidExtReq,Q,mdotBetweenTanks,mdotVent,QInclFromEnv] = findControlPointForConditions(obj,mdotFillRate,internalTempChangeRate,externalTempChangeRate,externalFluidTempChangeRate,isMaxingIntTempChangeRate)
           if ~exist('isMaxingIntTempChangeRate','var')
              isMaxingIntTempChangeRate = false; 
           end
           [mdotBetweenTanks,mdotVent,TFluidExtReq,Q,QInclFromEnv] = obj.findPointForConditions(mdotFillRate,internalTempChangeRate,externalTempChangeRate,externalFluidTempChangeRate,isMaxingIntTempChangeRate);
           if(mdotVent < 0 || obj.internalTank.vapourPressure < FillingSystem.ATM_PRESSURE+500e3 || isMaxingIntTempChangeRate) %Requested internal tank temp change rate too high, just don't vent
               mdotVent = 0;
               mdotBetweenTanks = mdotFillRate;
           end
           if(mdotBetweenTanks == 0)
               fillValveOpenAmt = 0;
           else
                fillValveOpenAmt = obj.findFillValvePositionForFlow(mdotBetweenTanks);
           end
           if(mdotVent == 0)
              ventValveOpenAmt = 0; 
           else
              ventValveOpenAmt = obj.findVentValvePositionForFlow(mdotVent);
           end
       end
       
       function [fillValveOpenAmt,ventValveOpenAmt,TFluidExtReq,Q,mdotBetweenTanks,mdotVent,QInclFromEnv] = findControlEquilibriumPointForFillRate(obj,mdotFillRate)
           [mdotBetweenTanks,mdotVent,TFluidExtReq,Q,QInclFromEnv] = obj.findEquilibriumPointForFillRate(mdotFillRate);
           %Find fillValveOpenAmt that gives mdotBetweenTanks
           
%            disp("Between tanks: "+mdotBetweenTanks);
%            disp("Vent: "+mdotVent);
%            disp("T Ext fluid: "+TFluidExtReq);
%            disp("Q: "+Q);
           
           fillValveOpenAmt = obj.findFillValvePositionForFlow(mdotBetweenTanks);
           ventValveOpenAmt = obj.findVentValvePositionForFlow(mdotVent);
%             try
%            fillValveOpenAmt = betterfzero(@getMDotBetweenTanksErr,0.5,0,1,1e-7,1000,0.5*10^-16,true);
%             catch excep
%                 disp("Error finding fill valve open amt for mass flow "+mdotBetweenTanks);
%                 rethrow excep;
%             end
%             try
%                 ventValveOpenAmt = betterfzero(@getMDotVentErr,0.5,0,1,1e-7,1000,0.5*10^-16,true);
%             catch excep
%                 disp("Error finding vent valve open amt for mass flow "+mdotVent);
%                 rethrow excep;
%             end
%            
%            function mdotBetweenTanksErr = getMDotBetweenTanksErr(fillValveOpenAmt)
%                mdotBetweenTanksErr = mdotBetweenTanks - obj.getFlowBetweenTanksIfFillValveOpenAmtWasThis(fillValveOpenAmt);
%            end
%            
%            function mdotVentErr = getMDotVentErr(ventValveOpenAmt)
%                mdotVentErr = mdotVent - obj.getFlowOutIfVentValveOpenAmtWasThis(ventValveOpenAmt);
%            end
       end
       
       function [mdotBetweenTanks,mdotVent,TFluidExtReq,Q,QIncludingFromEnv] = findPointForConditions(obj,mdotFillRate,internalTempChangeRate,externalTempChangeRate,externalFluidTempChangeRate,isMaxingIntTempChangeRate)
           if ~exist('isMaxingIntTempChangeRate','var')
              isMaxingIntTempChangeRate = false; 
           end
           if isMaxingIntTempChangeRate
                mdotVent = 0;
           else
                mdotVent = obj.calcMDotOutThisInternalTankReqInclTempChangeRate(mdotFillRate,internalTempChangeRate);
           end
           mdotBetweenTanks = mdotVent + mdotFillRate;
           %Calculate Q
           QTankReq = obj.externalTankSystem.externalTank.findHeatRateInputForTempChangeRateWithLiquidDrainingRate(externalTempChangeRate,mdotBetweenTanks);
           TFluidExtReq = obj.externalTankSystem.calcFluidTempRequiredForHeatFluxToTank(QTankReq);
           [Q,QIncludingFromEnv] = obj.externalTankSystem.calcExternalHeatFluxReqForFluidTempChangeRate(TFluidExtReq,obj.ambientTemp,externalFluidTempChangeRate);
       end
       
       function [mdotBetweenTanks,mdotVent,TFluidExtReq,Q,QIncludingFromEnv] = findEquilibriumPointForFillRate(obj,mdotFillRate)
           mdotVent = obj.calcMDotOutThisInternalTankReq(mdotFillRate);
           mdotBetweenTanks = mdotVent + mdotFillRate;
           %Calculate Q
           QTankReq = obj.externalTankSystem.externalTank.findHeatRateInputToKeepTempForDrainingLiquidRate(mdotBetweenTanks);
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
       
       function DExtTdt = calcDExtTdtForDiffExtFluidT(obj,fillValveOpenAmt,TFluidExt,QExt)
           dt = 0.25;
           clone = copy(obj);
           T1 = obj.externalTank.temp;
           mdotFill = obj.pipeBetweenTanks.getMassFlow(obj.externalTank.vapourPressure,obj.internalTank.vapourPressure,fillValveOpenAmt);
           EFlowBetweenTanks = dt.* mdotFill .* FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.externalTank.temp,obj.externalTank.vapourPressure);
           clone.externalTankSystem.TFluid = TFluidExt;
           QExtTankToApply = dt .* obj.externalTankSystem.externalTankSurfaceArea .* obj.externalTankSystem.heatTransferCoeffTankWithFluid .* (TFluidExt - obj.externalTank.temp);
           clone.externalTank.changeNitrousMassEnergy(-mdotFill.*dt,-EFlowBetweenTanks+QExtTankToApply);
           DExtTdt = (clone.externalTank.temp-T1) ./ dt;
       end
       
       function advanceSystemBySmallTimeIncrem(obj,dt)
            mdotFill = obj.pipeBetweenTanks.getMassFlow(obj.externalTank.vapourPressure,obj.internalTank.vapourPressure,obj.fillValveOpenAmt);
            mdotVent = obj.ventPipe.getMassFlow(obj.internalTank.vapourPressure,FillingSystem.ATM_PRESSURE,obj.ventValveOpenAmt);
            EFlowBetweenTanks = dt.* mdotFill .* FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.externalTank.temp,obj.externalTank.vapourPressure);
            if(obj.internalTank.liquidHeight >= obj.internalTank.gasVentHoleHeight)
               hIntVenting = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(obj.internalTank.temp,obj.internalTank.vapourPressure);
            else
               hIntVenting = FluidType.NITROUS_GAS.getSpecificEnthalpy(obj.internalTank.temp,obj.internalTank.vapourPressure); 
            end
            EFlowOutOfRunTank = dt.* mdotVent .* hIntVenting;
            %obj.internalTank.changeNitrousMassEnergy(-mdotVent.*dt,-EFlowOutOfRunTank);
            %HEAT TRANSFERS
            QInternal = dt .* obj.calcCurrentHeatTransferRateAmbientToInternalTank();
            %obj.internalTank.addHeat(QInternal);
            
            dmInt = mdotFill.*dt-mdotVent.*dt;
            if(obj.internalTank.mTotalNitrous + dmInt < 0) %Don't allow negative mass
               dmInt = -obj.internalTank.mTotalNitrous;
            end
%             disp("Old temp: "+obj.internalTank.temp);
            obj.internalTank.changeNitrousMassEnergy(dmInt,EFlowBetweenTanks-EFlowOutOfRunTank+QInternal);
%             disp("vent: "+mdotVent+" new temp: "+obj.internalTank.temp);
            QExtTankToApply = obj.externalTankSystem.doHeatTransferForTimeStep(dt,obj.ambientTemp,obj.QExt);
            obj.externalTank.changeNitrousMassEnergy(-mdotFill.*dt,-EFlowBetweenTanks+QExtTankToApply);
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