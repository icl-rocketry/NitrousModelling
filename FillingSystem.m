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
               externalTankFluidHeatTransferCoeffWithSurroundings,mFluid)
            obj.ambientTemp = ambientTemp;
            obj.internalTank = GeometricNitrousAndAirTank(ambientTemp,initialInternalNitrousMass,0,internalTankHeight,internalTankCrossSectionA,internalVentHoleHeight); 
            obj.externalTank = GeometricNitrousAndAirTank(ambientTemp,initialExternalNitrousMass,0,externalTankHeight,externalTankCrossSectionA,externalTankHeight);
            obj.ventPipe = SaturatedPipeValvePipeFast('preBakedData/saturatedGasPipeValveFlowRates.mat');
            obj.pipeBetweenTanks = SaturatedPipeValvePipeFast('preBakedData/saturatedLiquidPipeValveFlowRates.mat');
            obj.internalTankSurfaceArea = internalTankSurfaceArea;
            obj.heatTransferCoeffInternalTank = heatTransferCoeffInternalTank;
            obj.externalTankSystem = ExternalTankSystem(obj.externalTank,externalTankSurfaceArea,...
                heatTransferCoeffExternalTankWithFluid,ambientTemp, SHCFluid, ...
               externalTankFluidSurfaceArea,externalTankFluidHeatTransferCoeffWithSurroundings,...
               mFluid);
       end
       
       function xDot = getSystemChangeRateVector(obj)
            x1 = obj.getStateVector();
            obj2 = copy(obj);
            dt = 0.01;
            obj2.advanceSystemBySmallTimeIncrem(dt);
            x2 = obj2.getStateVector();
            xDot = (x2-x1)./dt;
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
            QInternal = dt .* obj.internalTankSurfaceArea .* obj.heatTransferCoeffInternalTank .* (obj.ambientTemp - obj.internalTank.temp);
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
           obj2.TExternalFluid = x(6);
           obj2.fillValveOpenAmt = u(1);
           obj2.ventValveOpenAmt = u(2);
           obj2.QExt = u(3);
       end
       
       function x = getStateVector(obj)
           x = [obj.ambientTemp;
               obj.externalTank.temp;
               obj.externalTank.mTotalNitrous;
               obj.internalTank.temp;
               obj.internalTank.mTotalNitrous];
       end
       
       function u = getControlVector(obj)
           u = [obj.fillValveOpenAmt;
               obj.ventValveOpenAmt;
               obj.QExt];
       end
   end
end