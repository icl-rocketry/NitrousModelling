classdef FillingSystemTrajectoryGen
   properties
      fillingPreferredInternalTankTemp = -20+273.15;
      fillingPreferredExternalTankTemp = -18+273.15;
      finalTargetInternalTankTemp = 30+273.15;
      maxQExt = 2000; %Max 2000W heating ability
      minQExt = -100; %Max 100W cooling ability
      MAX_T_EXT_TANK = SaturatedNitrous.T_CRIT - 3;
      targetInternalTankEndMass = 8; %8Kg
      dtStep = 0.25; %sec per time increment
      maxPreferredFillingRate = 100e-3; %100g/sec
      maxQChangeRate = 50; %50W/sec max change rate of Q
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
       
       function grad = gradientForSmoothApproach(obj,setPointWant,currentValue,maxGradPositive,maxGradNegative,smoothTime,aggressivenessExponent)
           if ~exist('smoothTime','var')
              smoothTime = 5; 
           end
           if ~exist('aggressivenessExponent','var')
              aggressivenessExponent = 2; 
           end
           dVal = setPointWant-currentValue;
           if(dVal > 0)
              maxGrad = maxGradPositive; 
           else
              maxGrad = maxGradNegative;
           end
           %Smooth over what at max gradient would take 5 sec
           dValForSmoothing = maxGrad .* smoothTime; %dT ~ (dT/dt)*dt
           if(abs(dVal) > abs(dValForSmoothing))
               grad = maxGrad;
               return;
           end

           if (dVal > 0 && maxGrad < 0) || (dVal < 0 && maxGrad > 0)
               grad = maxGrad;
               return;
           end
           %Decrease gradient to zero as approach point
           SF = (abs(dVal) / abs(dValForSmoothing))^(1./aggressivenessExponent);
           grad = maxGrad * SF;
           if(grad > maxGradPositive)
               grad = maxGradPositive;
           elseif(grad < maxGradNegative)
                grad = maxGradNegative;
           end
       end
       
       function [dTIntdtMaxPositiveNow,dTIntdtMaxNegative] = getMaxTempChangeRatesCurrentlyPossible(obj,systemCurrent)
           [fillValveOpenAmt,dTIntdtMaxPositiveNow,dTextdt,dTextFdt] = systemCurrent.calcEffectOfFillingNoVentWithHeat(obj.maxPreferredFillingRate,0);
           dTIntdtMaxNegative = systemCurrent.calcEffectOfFullyVentingNoFill;
       end
       
       function u = determineControlInputToApply(obj,systemCurrent,xVector,uVector)
           %TODO Compare system current's state to desired end conditions and
           %imposed limits and decide what control input we want to apply
           %nominally (Respecting control imposed limits)
           currentTime = size(xVector,2)*obj.dtStep; %Eg. when xVector empty, time is zero
           prevI = size(xVector,2);
%            u = [1,0,0]; %Temporary for now
           TExt = systemCurrent.externalTank.temp;
           TInt = systemCurrent.internalTank.temp;
           Tfluid = systemCurrent.externalTankSystem.TFluid;
           [dTIntdtMaxPositiveNow,dTIntdtMaxNegative] = obj.getMaxTempChangeRatesCurrentlyPossible(systemCurrent);
           dTIntdtDesired = obj.gradientForSmoothApproach(obj.fillingPreferredInternalTankTemp,TInt,dTIntdtMaxPositiveNow,dTIntdtMaxNegative,5,0.8);
           isMaxingIntTempChangeRate = abs(dTIntdtDesired-dTIntdtMaxPositiveNow) < 0.01 || TInt < obj.fillingPreferredInternalTankTemp-1;
%            disp(dTIntdtDesired+" ("+dTIntdtMaxPositiveNow+", "+dTIntdtMaxNegative+") "+systemCurrent.internalTank.mTotalNitrous+" "+isMaxingIntTempChangeRate);
%            disp(systemCurrent.internalTank.temp);
%            disp(systemCurrent.internalTank.vapourPressure);
           [fillValveOpenAmt,ventValveOpenAmt,~,~,~,~,~] = systemCurrent.findControlPointForConditions(obj.maxPreferredFillingRate,dTIntdtDesired,0,0,isMaxingIntTempChangeRate);
           %Find max dTExtFdt achievable
           xdot1 = systemCurrent.calcEffectOfControlOnSystem(fillValveOpenAmt,ventValveOpenAmt,obj.maxQExt);
           dTfdtMax = xdot1(6);
           xdot2 = systemCurrent.calcEffectOfControlOnSystem(fillValveOpenAmt,ventValveOpenAmt,obj.minQExt);
           dTfdtMin = xdot2(6);
%            disp("Min dTf: "+dTfdtMin);
%            disp("Max dTf: "+dTfdtMax);
%            dTExtdt = xdot1(2); %dTExt/dt for current TExternalFluid
           highestTExternalFluidSoon = Tfluid + dTfdtMax .* 15; %The highest value of TExternalFluid possible 15 seconds later
           lowestTExternalFluidSoon = Tfluid + dTfdtMin .* 15; 
           %Pretend we have control even when we don't so don't end up
           %doing the opposite of what need to
           if(lowestTExternalFluidSoon > Tfluid)
               lowestTExternalFluidSoon = Tfluid-1;
           end
           if(highestTExternalFluidSoon < Tfluid)
               highestTExternalFluidSoon = Tfluid+1;
           end
%            disp("Lowest T external fluid: "+lowestTExternalFluidSoon);
           dTExtdtMax = systemCurrent.calcDExtTdtForDiffExtFluidT(fillValveOpenAmt,highestTExternalFluidSoon,0);
           dTExtdtMin = systemCurrent.calcDExtTdtForDiffExtFluidT(fillValveOpenAmt,lowestTExternalFluidSoon,0);
%            disp("Predicting dTExt/dt with fill valve pos: "+fillValveOpenAmt+" and vent pos: "+ventValveOpenAmt);
%            xdot3 = systemCurrent.calcEffectOfControlOnSystemWithDiffExtFluidT(fillValveOpenAmt,ventValveOpenAmt,obj.maxQExt,highestTExternalFluidSoon);
%            dTExtdtMax = xdot3(2);
%            xdot4 = systemCurrent.calcEffectOfControlOnSystemWithDiffExtFluidT(fillValveOpenAmt,ventValveOpenAmt,obj.minQExt,lowestTExternalFluidSoon);
%            dTExtdtMin = xdot4(2);
%            disp("dTExt/dt max: "+dTExtdtMax+" dTExt/dt min: "+dTExtdtMin+" cur TExt: "+TExt+" (aim for "+obj.fillingPreferredExternalTankTemp+")");
           dTExtDesired = obj.gradientForSmoothApproach(obj.fillingPreferredExternalTankTemp,TExt,dTExtdtMax,dTExtdtMin,12,0.9);
%            disp("dTExt des: "+dTExtDesired+", max: "+dTExtdtMax+" min: "+dTExtdtMin);
           %            if(prevI > 2)
%                 dTExtdtPrev = (xVector(2,prevI) - xVector(2,prevI-1)) / obj.dtStep;
%                 dTExt = obj.fillingPreferredExternalTankTemp - TExt;
%                 if((dTExt < 0 && dTExtdtPrev < 0) || (dTExt > 0 && dTExtdtPrev > 0))
%                    %If previous trajectory of ext tank temp is in correct direction
%                    
%                 end
%            end
           
           [fillValveOpenAmt,ventValveOpenAmt,TFluidExtReq,~,~,~,~] = systemCurrent.findControlPointForConditions(obj.maxPreferredFillingRate,dTIntdtDesired,dTExtDesired,0,isMaxingIntTempChangeRate);
%            disp("Actual fill valve pos: "+fillValveOpenAmt);
           drawnow;
           dTExtFDesired = obj.gradientForSmoothApproach(TFluidExtReq,Tfluid,dTfdtMax,dTfdtMin,2,0.9);
           [~,~,~,Q,~,~,~] = systemCurrent.findControlPointForConditions(obj.maxPreferredFillingRate,dTIntdtDesired,dTExtDesired,dTExtFDesired,isMaxingIntTempChangeRate);
%            disp("Want dTExt/dt: "+dTExtDesired+" which needs Tfluid: "+TFluidExtReq+" (now "+Tfluid+", tank T: "+TExt+") cur dTextF/dt wanted: "+dTExtFDesired+" Q: "+Q);
           if(Q > obj.maxQExt)
              Q = obj.maxQExt; 
           end
           if(Q < obj.minQExt)
              Q = obj.minQExt; 
           end
           if(currentTime > 0)
               %Impose a limit on change rate of Q
               dQ = Q - uVector(3,prevI);
               if(abs(dQ) > obj.maxQChangeRate*obj.dtStep) %TODO Make this limit on change rate not dependent on sim rate
                  dQ = obj.maxQChangeRate.*obj.dtStep .* (dQ/abs(dQ)); 
               end
               Q = uVector(3,prevI) + dQ;
           end
           u = [fillValveOpenAmt,ventValveOpenAmt,Q];
       end
       
       function [xVector,uVector] = propogateTrajectory(obj,systemCurrent,xVector,uVector)
           %TODO IF trajectory end point reached, then return
           currentTime = size(xVector,2)*obj.dtStep; %Eg. when xVector empty, time is zero
           disp(currentTime);
           drawnow;
           if(currentTime > 120)
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