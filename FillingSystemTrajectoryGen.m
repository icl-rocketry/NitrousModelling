%Very crude attempt to come up with a first estimate for what valves,
%etc... should be doing at each point in time to achieve filling
classdef FillingSystemTrajectoryGen
   properties
      fillingPreferredInternalTankTemp = -20+273.15;
      fillingPreferredExternalTankTemp = -18+273.15;
      warmingPreferredExternalTankTemp = 33+273.15;
      finalTargetInternalTankTemp = 30+273.15;
      maxQExt = 2000; %Max 2000W heating ability
      minQExt = -100; %Max 100W cooling ability
      MAX_T_EXT_TANK = SaturatedNitrous.T_CRIT - 3;
      targetInternalTankEndMass = 8; %8Kg
      dtStep = 0.25; %sec per time increment
      dtStepWarming = 10; %10 sec per time increment (Very coarse)
      maxPreferredFillingRate = 100e-3; %100g/sec
      maxQChangeRateNotWarming = 50; %50W/sec max change rate of Q
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
           paramStrut.isWarmingToFinalTemp = false; %Initialize variable as a structure
           [x,u] = obj.propogateTrajectory(systemCopy,[],[],paramStrut);
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
       
       function u = determineWarmingControlInputToApply(obj,systemCurrent,xVector,uVector,coarse)
           if(~exist('coarse','var'))
               coarse = false;
           end
           currentTime = size(xVector,2)*obj.dtStepWarming; %Eg. when xVector empty, time is zero
           prevI = size(xVector,2);
           targetExtTankTemp = obj.warmingPreferredExternalTankTemp;
           %targetIntTankTemp = obj.finalTargetInternalTankTemp;
           TExt = systemCurrent.externalTank.temp;
           Tfluid = systemCurrent.externalTankSystem.TFluid;
           ventValveOpenAmt = 0;
           %Amount for correct fill rate
           fillValveOpenAmt = systemCurrent.findFillValvePositionForFlow(obj.maxPreferredFillingRate);
           %[fillValveOpenAmt,~,~,~,~,~,~] = systemCurrent.findControlPointForConditions(obj.maxPreferredFillingRate,0,0,0,true);
           
           if(TExt < targetExtTankTemp-2 || coarse)
               u = [fillValveOpenAmt,ventValveOpenAmt,obj.maxQExt];
               return;
           end
           
           xdot1 = systemCurrent.calcEffectOfControlOnSystem(fillValveOpenAmt,ventValveOpenAmt,obj.maxQExt);
           dTfdtMax = xdot1(6);
           xdot2 = systemCurrent.calcEffectOfControlOnSystem(fillValveOpenAmt,ventValveOpenAmt,obj.minQExt);
           dTfdtMin = xdot2(6);
           highestTExternalFluidSoon = Tfluid + dTfdtMax .* 15; %The highest value of TExternalFluid possible 15 seconds later
           lowestTExternalFluidSoon = Tfluid + dTfdtMin .* 15; 
           dTExtdtMax = systemCurrent.calcDExtTdtForDiffExtFluidT(fillValveOpenAmt,highestTExternalFluidSoon,0);
           dTExtdtMin = systemCurrent.calcDExtTdtForDiffExtFluidT(fillValveOpenAmt,lowestTExternalFluidSoon,0);
           dTExtDesired = obj.gradientForSmoothApproach(targetExtTankTemp,TExt,dTExtdtMax,dTExtdtMin,12,0.9);
           [fillValveOpenAmt,ventValveOpenAmt,TFluidExtReq,~,~,~,~] = systemCurrent.findControlPointForConditions(obj.maxPreferredFillingRate,0,dTExtDesired,0,true);
           dTExtFDesired = obj.gradientForSmoothApproach(TFluidExtReq,Tfluid,dTfdtMax,dTfdtMin,2,0.9);
           [~,~,~,Q,~,~,~] = systemCurrent.findControlPointForConditions(obj.maxPreferredFillingRate,0,dTExtDesired,dTExtFDesired,true);
%            disp("Want dTExt/dt: "+dTExtDesired+" which needs Tfluid: "+TFluidExtReq+" (now "+Tfluid+", tank T: "+TExt+") cur dTextF/dt wanted: "+dTExtFDesired+" Q: "+Q);
           if(Q > obj.maxQExt)
              Q = obj.maxQExt; 
           end
           if(Q < obj.minQExt)
              Q = obj.minQExt; 
           end
%            if(currentTime > 0)
%                %Impose a limit on change rate of Q
%                dQ = Q - uVector(3,prevI);
%                if(abs(dQ) > obj.maxQChangeRateNotWarming*obj.dtStepWarming) %TODO Make this limit on change rate not dependent on sim rate
%                   dQ = obj.maxQChangeRateNotWarming.*obj.dtStepWarming .* (dQ/abs(dQ)); 
%                end
%                Q = uVector(3,prevI) + dQ;
%            end
           u = [fillValveOpenAmt,ventValveOpenAmt,Q];
       end
       
       function [u,paramStrut] = determineControlInputToApply(obj,systemCurrent,xVector,uVector,paramStrut)
           if(systemCurrent.ambientTemp > obj.finalTargetInternalTankTemp)
                finalTargetExternalTankTemp = obj.finalTargetInternalTankTemp + 2;
           else
                finalTargetExternalTankTemp = obj.warmingPreferredExternalTankTemp;
           end
           
           wantedExtTankTemp = obj.fillingPreferredExternalTankTemp;
           wantedInternalTankTemp = obj.fillingPreferredInternalTankTemp;
           wantedFillRate = obj.gradientForSmoothApproach(obj.targetInternalTankEndMass,systemCurrent.internalTank.mTotalNitrous,obj.maxPreferredFillingRate,-obj.maxPreferredFillingRate,2,0.9);
           fractionOfMaxFillRate = (wantedFillRate ./ obj.maxPreferredFillingRate);
           proportionOfFinalFillManouver = 1-fractionOfMaxFillRate;
           if(proportionOfFinalFillManouver > 1)
               proportionOfFinalFillManouver = 1;
           end
           
           %TODO Compare system current's state to desired end conditions and
           %imposed limits and decide what control input we want to apply
           %nominally (Respecting control imposed limits)
           if(~paramStrut.isWarmingToFinalTemp)
               approxDistFromNeedingToWarmFastMethod = obj.estimateWarmingFast(systemCurrent,obj.warmingPreferredExternalTankTemp);
               disp("Approx warm distance thing: "+approxDistFromNeedingToWarmFastMethod);
               if(approxDistFromNeedingToWarmFastMethod < 100e3)
                   [distFromNeedingToWarm,~,~] = obj.propogateWarming(copy(systemCurrent),[],[]); %Most accurate, but v costly to work out
                   disp("Dist from needing to warm: "+distFromNeedingToWarm);
                   if(distFromNeedingToWarm < 0.2)
                       paramStrut.isWarmingToFinalTemp = true;
                   end
               end
           end
           
           if(paramStrut.isWarmingToFinalTemp)
               disp("SHOULD WARM");
               %TODO Smooth out final approach to end point
               wantedExtTankTemp = (1-proportionOfFinalFillManouver).*obj.warmingPreferredExternalTankTemp + proportionOfFinalFillManouver.*finalTargetExternalTankTemp;
               wantedInternalTankTemp = obj.finalTargetInternalTankTemp;
           end
           
           
           currentTime = size(xVector,2)*obj.dtStep; %Eg. when xVector empty, time is zero
           prevI = size(xVector,2);
%            u = [1,0,0]; %Temporary for now
           TExt = systemCurrent.externalTank.temp;
           TInt = systemCurrent.internalTank.temp;
           Tfluid = systemCurrent.externalTankSystem.TFluid;
           [dTIntdtMaxPositiveNow,dTIntdtMaxNegative] = obj.getMaxTempChangeRatesCurrentlyPossible(systemCurrent);
           dTIntdtDesired = obj.gradientForSmoothApproach(wantedInternalTankTemp,TInt,dTIntdtMaxPositiveNow,dTIntdtMaxNegative,5,0.8);
           isMaxingIntTempChangeRate = abs(dTIntdtDesired-dTIntdtMaxPositiveNow) < 0.01 || TInt < obj.fillingPreferredInternalTankTemp-1;
%            disp(dTIntdtDesired+" ("+dTIntdtMaxPositiveNow+", "+dTIntdtMaxNegative+") "+systemCurrent.internalTank.mTotalNitrous+" "+isMaxingIntTempChangeRate);
%            disp(systemCurrent.internalTank.temp);
%            disp(systemCurrent.internalTank.vapourPressure);
           [fillValveOpenAmt,ventValveOpenAmt,~,~,~,~,~] = systemCurrent.findControlPointForConditions(wantedFillRate,dTIntdtDesired,0,0,isMaxingIntTempChangeRate);
           %Find max dTExtFdt achievable
           xdot1 = systemCurrent.calcEffectOfControlOnSystem(fillValveOpenAmt,ventValveOpenAmt,obj.maxQExt);
           dTfdtMax = xdot1(6);
           xdot2 = systemCurrent.calcEffectOfControlOnSystem(fillValveOpenAmt,ventValveOpenAmt,obj.minQExt);
           dTfdtMin = xdot2(6);
%            disp("Min dTf: "+dTfdtMin);
%            disp("Max dTf: "+dTfdtMax);
%            dTExtdt = xdot1(2); %dTExt/dt for current TExternalFluid
           highestTExternalFluidSoon = Tfluid + dTfdtMax .* 30; %The highest value of TExternalFluid possible 15 seconds later
           lowestTExternalFluidSoon = Tfluid + dTfdtMin .* 30; 
           %Pretend we have control even when we don't so don't end up
           %doing the opposite of what need to
           if(lowestTExternalFluidSoon > Tfluid)
               lowestTExternalFluidSoon = Tfluid-2;
           end
           if(highestTExternalFluidSoon < Tfluid)
               highestTExternalFluidSoon = Tfluid+2;
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
           dTExtDesired = obj.gradientForSmoothApproach(wantedExtTankTemp,TExt,dTExtdtMax,dTExtdtMin,12,0.9);
%            disp("dTExt des: "+dTExtDesired+", max: "+dTExtdtMax+" min: "+dTExtdtMin);
           %            if(prevI > 2)
%                 dTExtdtPrev = (xVector(2,prevI) - xVector(2,prevI-1)) / obj.dtStep;
%                 dTExt = obj.fillingPreferredExternalTankTemp - TExt;
%                 if((dTExt < 0 && dTExtdtPrev < 0) || (dTExt > 0 && dTExtdtPrev > 0))
%                    %If previous trajectory of ext tank temp is in correct direction
%                    
%                 end
%            end
           
           [fillValveOpenAmt,ventValveOpenAmt,TFluidExtReq,~,~,~,~] = systemCurrent.findControlPointForConditions(wantedFillRate,dTIntdtDesired,dTExtDesired,0,isMaxingIntTempChangeRate);
%            disp("Actual fill valve pos: "+fillValveOpenAmt);
           dTExtFDesired = obj.gradientForSmoothApproach(TFluidExtReq,Tfluid,dTfdtMax,dTfdtMin,2,0.9);
           [~,~,~,Q,~,~,~] = systemCurrent.findControlPointForConditions(wantedFillRate,dTIntdtDesired,dTExtDesired,dTExtFDesired,isMaxingIntTempChangeRate);
%            disp("Want dTExt/dt: "+dTExtDesired+" which needs Tfluid: "+TFluidExtReq+" (now "+Tfluid+", tank T: "+TExt+") cur dTextF/dt wanted: "+dTExtFDesired+" Q: "+Q);
           if(Q > obj.maxQExt)
              Q = obj.maxQExt; 
           end
           if(Q < obj.minQExt)
              Q = obj.minQExt; 
           end
           if(currentTime > 0 && ~paramStrut.isWarmingToFinalTemp)
               %Impose a limit on change rate of Q
               dQ = Q - uVector(3,prevI);
               if(abs(dQ) > obj.maxQChangeRateNotWarming*obj.dtStep) %TODO Make this limit on change rate not dependent on sim rate
                  dQ = obj.maxQChangeRateNotWarming.*obj.dtStep .* (dQ/abs(dQ)); 
               end
               Q = uVector(3,prevI) + dQ;
               if(Q > obj.maxQExt)
                   Q = obj.maxQExt;
               end
               if(Q < obj.minQExt)
                   Q = obj.minQExt;
               end
           end
%            disp("Final Q to apply: "+Q);
%            drawnow;
           u = [fillValveOpenAmt,ventValveOpenAmt,Q];
       end
       
       function distFromNeedingToWarm = estimateWarmingFast(obj,systemCurrent,TExtTankWarming)
           %Estimate temp difference between warming now and final target
           %temp
           mdotBetween = obj.maxPreferredFillingRate;
           dTExtTank = TExtTankWarming - systemCurrent.externalTank.temp;
           
           ERate = @(T) mdotBetween .* NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,T,FluidProperty.VAPOR_QUALITY,0);
           QToExtFluid = @(T) systemCurrent.externalTankSystem.fluidSurfaceAreaWithSurroundings .* systemCurrent.externalTankSystem.fluidHeatTransferCoeffWithSurroundings .* (systemCurrent.ambientTemp - T);
           %Pick a temp to use for the representative change rate of tank
           %temp to use to estimate
           avgExtTankWarmingTemp = 0.25.*(3.*systemCurrent.externalTank.temp+TExtTankWarming);
           ERateLoseFromExtTank = ERate(avgExtTankWarmingTemp); %Approximately
           QRateExtTank = obj.maxQExt + QToExtFluid(avgExtTankWarmingTemp) - ERateLoseFromExtTank; %VERY approximately
           extTankCopy = copy(systemCurrent.externalTank);
           extTankCopy.forceSetTemp(avgExtTankWarmingTemp);
           extTankCopy.changeNitrousMassEnergy(-mdotBetween.*0.2,QRateExtTank.*0.2);
           approxdExtTdt = (extTankCopy.temp-avgExtTankWarmingTemp) ./ 0.2;
%            disp("Average ext tank temp while heating it up: "+(avgExtTankWarmingTemp-273.15));
%            disp("mdot between tanks: "+mdotBetween);
           %            C = 0.5.*(NitrousFluid.getSaturatedHeatCapacity(systemCurrent.externalTank.temp)+NitrousFluid.getSaturatedHeatCapacity(TExtTankWarming)); %VERY approx
           %Q = mcdT
%            approxdExtTdt = QRateExtTank ./ (systemCurrent.externalTank.mTotalNitrous.*C);
           timeToBringToTExtTank = abs(dTExtTank ./ approxdExtTdt);
%            disp("Q Ext Tank from env through fluid: "+(QRateExtTank)+" Approx dExtTdt: "+approxdExtTdt+" time to bring to TExtTank: "+timeToBringToTExtTank);
           
           mCur = systemCurrent.internalTank.mTotalNitrous;
           dm = obj.targetInternalTankEndMass - mCur;
           t = dm ./ mdotBetween;
           if(timeToBringToTExtTank > t)
               distFromNeedingToWarm = 0; %TODO: Really ought to change the fill rate (Reduce it)
               return;
           end
           ERate2 = ERate(TExtTankWarming);
           ERate1 = 0.5.*(ERate2 + ERate(systemCurrent.externalTank.temp));
%            disp("Total time to fill: "+t);
%            disp("Time heating ext tank: "+timeToBringToTExtTank);
           dE = (t-timeToBringToTExtTank).*ERate1 + timeToBringToTExtTank.*ERate2;
           QRateInternalTank = @(T) systemCurrent.internalTankSurfaceArea .* systemCurrent.heatTransferCoeffInternalTank .* (systemCurrent.ambientTemp - T);
           %Assume average Q is just mean of start and end
           QIntTank = 0.5.*(QRateInternalTank(systemCurrent.internalTank.temp) + QRateInternalTank(obj.finalTargetInternalTankTemp));
           %            dQ = integral(QRateInternalTank,systemCurrent.internalTank.temp,obj.finalTargetInternalTankTemp);
           dE = dE + QIntTank.*t;
           tankCopy = copy(systemCurrent.internalTank);
           tankCopy.forceSetNitrousMass(obj.targetInternalTankEndMass);
           tankCopy.forceSetTemp(obj.finalTargetInternalTankTemp);
           dEReq = tankCopy.getInternalEnergy() - systemCurrent.internalTank.getInternalEnergy();
           distFromNeedingToWarm = dE - dEReq;
       end
       
       function [distFromNeedingToWarm,xVector,uVector] = propogateWarming(obj,systemCurrent,xVector,uVector)
           TInt = systemCurrent.internalTank.temp;
           if(systemCurrent.internalTank.mTotalNitrous >= obj.targetInternalTankEndMass)
              distFromNeedingToWarm = 0; %Will reach target end mass before warms up enough
              return;
           end
           if(TInt >= obj.finalTargetInternalTankTemp) %Tank has warmed to max
               %Distance from needing to warm is difference in mass between
               %final tank and target end mass. +ve distance is before need
               %to warm
               distFromNeedingToWarm = obj.targetInternalTankEndMass - systemCurrent.internalTank.mTotalNitrous;
               return;
           end
           u = obj.determineWarmingControlInputToApply(systemCurrent,xVector,uVector,true);
          
           systemCurrent.fillValveOpenAmt = u(1);
           systemCurrent.ventValveOpenAmt = u(2);
           systemCurrent.QExt = u(3);
%            disp((TInt-273.15)+" "+(systemCurrent.externalTank.temp-273.15)+" "+u(3)+" "+systemCurrent.externalTank.vapourPressure+" "+systemCurrent.externalTank.isSaturated+" "+systemCurrent.externalTank.mLiquid+" "+systemCurrent.externalTank.mVapour+" "+systemCurrent.internalTank.mTotalNitrous);
%            drawnow;
           %Update current system state into vector for this time step
           xVector(:,size(xVector,2)+1) = systemCurrent.getStateVector();
           uVector(:,size(uVector,2)+1) = systemCurrent.getControlVector();
           
           %Advance the system to the next time step and recursively call
           %ourself again
           systemCurrent.advanceSystemBySmallTimeIncrem(obj.dtStepWarming);
           
           [distFromNeedingToWarm,xVector,uVector] = obj.propogateWarming(systemCurrent,xVector,uVector);
       end
       
       function [xVector,uVector] = propogateTrajectory(obj,systemCurrent,xVector,uVector,paramStrut)
           %TODO IF trajectory end point reached, then return
           currentTime = size(xVector,2)*obj.dtStep; %Eg. when xVector empty, time is zero
           disp(currentTime);
           drawnow;
           if(currentTime > 340)
              return; %TODO REMOVE THIS TEMPORARY THING 
           end
           
           %Compare system current's state to desired end conditions and
           %imposed limits and decide what control input we want to apply
           %nominally (Respecting control imposed limits)
           [u,paramStrut] = obj.determineControlInputToApply(systemCurrent,xVector,uVector,paramStrut);
           systemCurrent.fillValveOpenAmt = u(1);
           systemCurrent.ventValveOpenAmt = u(2);
           systemCurrent.QExt = u(3);
           
           %Update current system state into vector for this time step
           xVector(:,size(xVector,2)+1) = systemCurrent.getStateVector();
           uVector(:,size(uVector,2)+1) = systemCurrent.getControlVector();
           
           %Advance the system to the next time step and recursively call
           %ourself again
           systemCurrent.advanceSystemBySmallTimeIncrem(obj.dtStep);
           [xVector,uVector] = obj.propogateTrajectory(systemCurrent,xVector,uVector,paramStrut);
       end
   end
end