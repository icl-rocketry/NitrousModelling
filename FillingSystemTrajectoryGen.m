classdef FillingSystemTrajectoryGen
   properties
      fillingPreferredInternalTankTemp = -20+273.15;
      fillingPreferredExternalTankTemp = -18+273.15;
      finalTargetInternalTankTemp = 30+273.15;
      maxQExt = 2000; %Max 2000W heating ability
      minQExt = -500; %Max 500W cooling ability
      MAX_T_EXT_TANK = SaturatedNitrous.T_CRIT - 3;
      targetInternalTankEndMass = 8; %8Kg
      dtStep = 0.1; %0.1 sec per time increment
      maxPreferredFillingRate = 100e-3; %100g/sec
   end
   
   methods
       function obj = FillingSystemTrajectoryGen(fillingPreferredInternalTankTemp,...
               fillingPreferredExternalTankTemp,finalTargetInternalTankTemp,...
               maxQExt,minQExt,targetInternalTankEndMass,maxPreferredFillingRate)
           obj.fillingPreferredInternalTankTemp = fillingPreferredInternalTankTemp;
           obj.fillingPreferredExternalTankTemp = fillingPreferredExternalTankTemp;
           obj.finalTargetInternalTankTemp = finalTargetInternalTankTemp;
           obj.maxQExt = maxQExt;
           obj.minQExt = minQExt;
           obj.targetInternalTankEndMass = targetInternalTankEndMass;
           obj.maxPreferredFillingRate = maxPreferredFillingRate;
       end
       
       %Returns instance of FillingSystemTrajectory
       function traj = genTrajectory(obj,system)
           systemCopy = copy(system);
           [x,u] = obj.propogateTrajectory(systemCopy,[],[]);
           traj = FillingSystemTrajectory(obj.dtStep,x,u);
       end
       
       function u = determineControlInputToApply(obj,systemCurrent,xVector,uVector)
           %TODO Compare system current's state to desired end conditions and
           %imposed limits and decide what control input we want to apply
           %nominally (Respecting control imposed limits)
           currentTime = size(xVector,2)*obj.dtStep; %Eg. when xVector empty, time is zero
           u = [1,0,0]; %Temporary for now
       end
       
       function [xVector,uVector] = propogateTrajectory(obj,systemCurrent,xVector,uVector)
           %TODO IF trajectory end point reached, then return
           currentTime = size(xVector,2)*obj.dtStep; %Eg. when xVector empty, time is zero
           disp(currentTime);
           drawnow;
           if(currentTime > 60)
              return; %TODO REMOVE THIS TEMPORARY THING 
           end
           
           %Compare system current's state to desired end conditions and
           %imposed limits and decide what control input we want to apply
           %nominally (Respecting control imposed limits)
           u = obj.determineControlInputToApply(systemCurrent,xVector,uVector);
           systemCurrent.fillValveOpenAmt = u(1);
           systemCurrent.ventValveOpenAmt = u(2);
           systemCurrent.QExt = u(3);
           
           %Update current system state into vector for this time step
           xVector(:,size(xVector,2)+1) = systemCurrent.getStateVector();
           uVector(:,size(uVector,2)+1) = systemCurrent.getControlVector();
           
           %Advance the system to the next time step and recursively call
           %ourself again
           systemCurrent.advanceSystemBySmallTimeIncrem(obj.dtStep);
           [xVector,uVector] = obj.propogateTrajectory(systemCurrent,xVector,uVector);
       end
   end
end