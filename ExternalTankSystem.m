%Represents the external tank submerged in a fluid whose temperature is
%being regulated by heat input/extraction
classdef ExternalTankSystem < matlab.mixin.Copyable
   properties
       externalTank;
       externalTankSurfaceArea;
       heatTransferCoeffTankWithFluid;
       TFluid;
       SHCFluid = 4200; %Specific heat capacity of fluid
       mFluid;
       fluidSurfaceAreaWithSurroundings;
       fluidHeatTransferCoeffWithSurroundings;
   end
   methods(Access = protected)
       % Override copyElement method:
       function cpObj = copyElement(obj)
           % Make a shallow copy of all primitive properties
           cpObj = copyElement@matlab.mixin.Copyable(obj);
           % Do a deep copy of handle properties
           cpObj.externalTank = copy(obj.externalTank);
       end
   end
   methods
       function obj = ExternalTankSystem(externalTank,externalTankSurfaceArea, ...
               heatTransferCoeffTankWithFluid,TFluidInitial, SHCFluid, ...
               fluidSurfaceAreaWithSurroundings,fluidHeatTransferCoeffWithSurroundings,...
               mFluid)
          obj.externalTank = externalTank; 
          obj.externalTankSurfaceArea = externalTankSurfaceArea;
          obj.heatTransferCoeffTankWithFluid = heatTransferCoeffTankWithFluid;
          obj.TFluid = TFluidInitial;
          obj.SHCFluid = SHCFluid;
          obj.fluidSurfaceAreaWithSurroundings = fluidSurfaceAreaWithSurroundings;
          obj.fluidHeatTransferCoeffWithSurroundings = fluidHeatTransferCoeffWithSurroundings;
          obj.mFluid = mFluid;
       end
       
       function TFluid = calcFluidTempRequiredForHeatFluxToTank(obj,Q)
           TFluid = obj.externalTank.temp + (Q/(obj.externalTankSurfaceArea .* obj.heatTransferCoeffTankWithFluid));
       end
       
       function [QExt,QTotalReq] = calcExternalHeatFluxReqForFluidTempChangeRate(obj,TFluidReq,TAmbient,dTFdt)
           QToTank = obj.externalTankSurfaceArea .* obj.heatTransferCoeffTankWithFluid .* (TFluidReq - obj.externalTank.temp);
           QAmbient = obj.fluidSurfaceAreaWithSurroundings .* obj.fluidHeatTransferCoeffWithSurroundings .* (TAmbient - TFluidReq);
           QChangeT = obj.mFluid.*obj.SHCFluid.*(dTFdt); %Q = mc(dT)
           %QReq = QExt + QAmbient
           QExt = QToTank - QAmbient + QChangeT;
           QTotalReq = QToTank + QChangeT;
       end
       
       function [QExt,QTotalReq] = calcExternalHeatFluxReqToKeepFluidAtTemp(obj,TFluidReq,TAmbient)
           QReq = obj.externalTankSurfaceArea .* obj.heatTransferCoeffTankWithFluid .* (TFluidReq - obj.externalTank.temp);
           QAmbient = obj.fluidSurfaceAreaWithSurroundings .* obj.fluidHeatTransferCoeffWithSurroundings .* (TAmbient - TFluidReq);
           %QReq = QExt + QAmbient
           QExt = QReq - QAmbient;
           QTotalReq = QReq;
       end
       
       %Note Q is energy (J) not power (J/s) here
       function applyHeatToFluid(obj,Q)
          dT = Q./(obj.mFluid.*obj.SHCFluid); %Q = mc(Dt)
          obj.TFluid = obj.TFluid + dT;
       end
       
       %Note Q is power (J/s) not energy (J) here
       function doHeatTransferForTimeStep(obj,dt,TAmbient,QInputToFluid)
           Q1 = dt .* obj.fluidSurfaceAreaWithSurroundings .* obj.fluidHeatTransferCoeffWithSurroundings .* (TAmbient - obj.TFluid);
           Q2 = dt .* QInputToFluid;
           obj.applyHeatToFluid(Q1+Q2);
           QTank = dt .* obj.externalTankSurfaceArea .* obj.heatTransferCoeffTankWithFluid .* (obj.TFluid - obj.externalTank.temp);
           obj.externalTank.addHeat(QTank);
           obj.applyHeatToFluid(-QTank);
       end
   end
end