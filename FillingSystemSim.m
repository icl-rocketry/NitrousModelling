clear
clc
close all

ambientTemp = 0+273.15;
internalTankHeight = 0.8;
internalTankDiameter = 150e-3;
externalTankDiameter = 400e-3;
internalTankCrossSectionA = 0.25*pi*(internalTankDiameter).^2;
internalVentHoleHeight = 0.8.*internalTankHeight;
initialInternalNitrousMass = 4;
initialExternalNitrousMass = 50;
externalTankHeight = 1;
externalTankCrossSectionA = 0.25*pi*(externalTankDiameter).^2;
internalTankSurfaceArea = pi*internalTankDiameter*internalTankHeight+2*internalTankCrossSectionA;
internalTankWallThickness = 2e-3;
externalTankWallThickness = 4e-3; %Assumed
% %Overall heat transfer coeff from just Stainless 316 thermal conductivity, is an oversimplification
% %Thermal conductivity is W/mK, to make into overall coefficient (W/m^2K)
% %need to divide by thickness
% heatTransferCoeffInternalTank = 13.5 / internalTankWallThickness;

%Heat transfer coeff stainless steel in series with some insulation
%material. BIG oversimplification
insulationThickness = 0e-3;%5e-3;
airConvectionHeatTransferCoeff = 10.45;%10.45; %TODO
nitrousConvectionHeatTransferCoeff = 44154; %TODO

%(1/nitrousConvectionHeatTransferCoeff)
heatTransferCoeffInternalTank = 1/((1/airConvectionHeatTransferCoeff) + 0 + (insulationThickness/0.03) + (internalTankWallThickness/13.5) );

externalTankSurfaceArea = pi*externalTankDiameter*externalTankHeight+2*externalTankCrossSectionA;
%Overall heat transfer coeff jusst from Stainless 316 thermal conductivity,
%is an oversimplifcation
heatTransferCoeffExternalTankWithFluid = 13.5 /  externalTankWallThickness; %Stainless 316 assumed
SHCFluid = 3000; %Very approx
%Vague guess
externalTankFluidSurfaceArea = pi*(externalTankDiameter+2e-2)*externalTankHeight+2*0.25*pi*(externalTankDiameter+2e-2)^2;
%Bad estimate, but assuming simply 15mm of substance with thermal
%conductivity 0.2 and only conduction
externalTankFluidHeatTransferCoeffWithSurroundings = 0.2/15e-3;
%Simplistic guess
extFluidVol = externalTankHeight*(0.25*pi*(externalTankDiameter+2e-2)^2-0.25*pi*(externalTankDiameter)^2);
mFluid = 997*extFluidVol; %Water density of 997 assumed

initialInternalTankTemp = 27+273.15;
initialExternalTankTemp = 34+273.15;
initialExtFluidTemp = 34+273.15;
% initialExtFluidTemp = -17.5746+273.15;
% initialInternalTankTemp = -20+273.15;
% initialExternalTankTemp = 34+273.15;
% initialExtFluidTemp = 34+273.15;

fillingSystem = FillingSystem(ambientTemp,internalTankHeight,internalTankCrossSectionA...
               ,internalVentHoleHeight,initialInternalNitrousMass,initialExternalNitrousMass,...
               externalTankHeight,externalTankCrossSectionA,internalTankSurfaceArea,...
               heatTransferCoeffInternalTank,externalTankSurfaceArea,...
               heatTransferCoeffExternalTankWithFluid,SHCFluid,externalTankFluidSurfaceArea,...
               externalTankFluidHeatTransferCoeffWithSurroundings,mFluid,initialInternalTankTemp,...
               initialExternalTankTemp,initialExtFluidTemp);

disp("Ext tank initial state: <Is saturated> <Temp> <MVapour> <MLiquid>");
disp(fillingSystem.externalTank.isSaturated+" "+(fillingSystem.externalTank.temp-273.15)+" "+fillingSystem.externalTank.mVapour+" "+fillingSystem.externalTank.mLiquid);
           
mdotFillRate = 40e-3;
internalTempChangeRate = 0 / 60.0; %C/min -> C/sec
externalTempChangeRate = 0 / 60.0; %C/sec
externalFluidTempChangeRate = 0 / 60.0; %C/sec
[fillValveOpenAmt,ventValveOpenAmt,TFluidExtReq,Q,mdotBetweenTanks,mdotVent,QInclFromEnv] = fillingSystem.findControlPointForConditions(mdotFillRate,internalTempChangeRate,externalTempChangeRate,externalFluidTempChangeRate);
disp("Heat transfer overall coefficient internal tank with surroundings: "+heatTransferCoeffInternalTank+" W/m^2K");
disp("Heat transfer overall coefficient external fluid with surroundings: "+externalTankFluidHeatTransferCoeffWithSurroundings+" W/m^2K");
disp("Heat transfer overall coefficient external tank with fluid: "+heatTransferCoeffExternalTankWithFluid+" W/m^2K");
disp("Ambient temp: "+(fillingSystem.ambientTemp-273.15)+" C");
disp("Internal tank temp: "+(fillingSystem.internalTank.temp-273.15)+" C");
disp("Internal tank temp change rate: "+(internalTempChangeRate*60)+" C/min");
disp("Internal tank pressure: "+fillingSystem.internalTank.vapourPressure);
disp("External tank temp: "+(fillingSystem.externalTank.temp-273.15)+" C");
disp("External tank temp change rate: "+(externalTempChangeRate*60)+" C/min");
disp("Fill rate of "+(mdotFillRate*1000)+" g/sec");
disp("Fill valve position: "+(fillValveOpenAmt*100)+"% open");
disp("Mdot due to fill valve: "+(mdotBetweenTanks*1000)+" g/sec");
disp("Vent valve position: "+(ventValveOpenAmt*100)+"% open");
disp("Mdot due to vent valve: "+(mdotVent*1000)+" g/sec");
disp("T Ext fluid: "+(TFluidExtReq-273.15)+" C");
disp("T Ext fluid change rate: "+(externalFluidTempChangeRate*60)+" C/min");
disp("Required Extra Q flux into Ext fluid (Eg. on top of heat transfer from ambient environment): "+Q+" W");
disp("Required total Q flux into Ext fluid: "+QInclFromEnv+" W");
nitrousEfficiency = (mdotFillRate*100) / mdotBetweenTanks;
disp("Nitrous efficiency: "+nitrousEfficiency+"%");

disp("----------------------------");
[fillValveOpenAmt2,dTintdt2,dTextdt2,dTextFdt2] = fillingSystem.calcEffectOfFillingNoVentWithHeat(mdotFillRate,0);
disp("If close vent valve and stop applying any heat:");
disp("To maintain fill rate of "+(mdotFillRate*1000)+" g/sec need fill valve "+(fillValveOpenAmt2*100)+"% open");
disp("Internal tank temp change rate: "+(dTintdt2*60)+" C/min");
disp("External tank temp change rate: "+(dTextdt2*60)+" C/min");
disp("T Ext fluid used: "+(fillingSystem.externalTankSystem.TFluid-273.15)+" C");
disp("T Ext fluid change rate: "+(dTextFdt2*60)+" C/min");

disp("----------------------------");
mdotFillRate = 40e-3;
[fillValveOpenAmt3,dTintdt3,TFluidExtReq3,Q3,QInclFromEnv3] = fillingSystem.calcFillingPointNoVentExternalTankEquilibrium(mdotFillRate);
disp("If close vent valve and want to stabilize external tank system:");
disp("To maintain fill rate of "+(mdotFillRate*1000)+" g/sec need fill valve "+(fillValveOpenAmt3*100)+"% open");
disp("External tank fluid temp required: "+(TFluidExtReq3-273.15)+" C");
disp("Internal tank temp change rate: "+(dTintdt3*60)+" C/min");
disp("Required Extra Q flux into Ext fluid (Eg. on top of heat transfer from ambient environment): "+Q3+" W");
disp("Required total Q flux into Ext fluid: "+QInclFromEnv3+" W");

disp("----------------------------");
fillingPreferredInternalTankTemp = -20+273.15;
fillingPreferredExternalTankTemp = 34+273.15;
finalTargetInternalTankTemp = 27+273.15;
maxQExt = 5000; %5000W max heating capacity
minQExt = -500; %Max 500W cooling ability
targetInternalTankEndMass = 5; %8Kg
maxPreferredFillingRate = 50e-3; %40g/sec
trajectoryGen = FillingSystemTrajectoryGen(fillingPreferredInternalTankTemp,...
               fillingPreferredExternalTankTemp,finalTargetInternalTankTemp,...
               maxQExt,minQExt,targetInternalTankEndMass,maxPreferredFillingRate);
           
fillingSystem.internalTank.forceSetNitrousMass(0);
fillingSystem.internalTank.forceSetTemp(SaturatedNitrous.getSaturationTemperature(101325));
           
tic;
trajectory = trajectoryGen.genTrajectory(fillingSystem);
toc;
vis = SystemVisualization();
while(true)
    t = 0;
    while(t<=340)
        t = t+0.1;
        [xSim,uSim] = trajectory.getXAndUForTime(t);
        vis.draw(xSim,uSim,t);
        pause(0.09);
    end
    pause(1);
end