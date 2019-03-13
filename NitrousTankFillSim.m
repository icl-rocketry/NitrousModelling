%Script to simulate filling of the nitrous tank
%Depends on all files with Nitrous in the name
clear
clc
CELSIUS_TO_K = 273.15; %Use a constant to help with code readability

ambientTemp = 28+CELSIUS_TO_K; %The ambient temp in K, also starting temp of internal tank

tankRadius = 0.065; %Radius of tank internal geometry in metres
tankHeight = 0.8;%1.35;%1.15; %Height of tank internal geometry in metres
global holeInInternalTankHeight;
holeInInternalTankHeight = 0.8*tankHeight;
tankCrossSectionA = pi * tankRadius^2;
internalTankSurfaceArea = 2*tankCrossSectionA + pi*2*tankRadius*tankHeight;
mAir = 1.225 * tankHeight * tankCrossSectionA; %Assume density of 1.225 for air starting in tank

%tank = GeometricNitrousAndAirTank(tankTemp,6,0,tankHeight,tankCrossSectionA);
runTank = GeometricNitrousAndAirTank(20+CELSIUS_TO_K,2,0,tankHeight,tankCrossSectionA,holeInInternalTankHeight);

% if true
%     
%    return; 
% end

global fillSimRunTank;
fillSimRunTank = runTank;

extTankRadius = 0.065;
extTankHeight = 3; %Height of tank external geometry in metres
extTankTemp = 28+CELSIUS_TO_K; %Temp to maintain external tank at
extTankCrossSectionA = pi * extTankRadius^2;
mNitrousExt = 10; %Not sure if physical

extTank = GeometricNitrousAndAirTank(extTankTemp,mNitrousExt,0,extTankHeight,extTankCrossSectionA,0);
global fillSimExtTank;
fillSimExtTank = extTank;

pipeBetweenTanksOutletCrossSectionA = pi * 0.25 * (0.5*4.8*10^-3)^2; %1/2 inch diameter
pipeBetweenTanksThinnestCrossSectionA = pi * 0.25 * (0.5*4.8*10^-3)^2; %1/2 inch diameter
pipeBetweenTanks = NitrousFluidPipe(@extTankConditions,@intTankInletPressure,pipeBetweenTanksOutletCrossSectionA,pipeBetweenTanksThinnestCrossSectionA);

ventPipeCrossSectionA = pi * 0.25 * (4.8*10^-3)^2; %1/4 inch diameter
%ventPipeThinnestCrossSectionA = ventPipeCrossSectionA;
ventPipeThinnestCrossSectionA = pi * 0.25 * (3*10^-3)^2; %0.3mm diameter
ventPipe = NitrousFluidPipe(@intTankConditions,@atmosphericPressure,ventPipeCrossSectionA,ventPipeThinnestCrossSectionA);

heatTransferCoeffInternal = 205; %Conductivity of aluminium, will be a large overestimate

fillingModel = NitrousTankFilling(runTank,extTank,pipeBetweenTanks,ventPipe,holeInInternalTankHeight,ambientTemp,heatTransferCoeffInternal,internalTankSurfaceArea,extTankTemp);

if false
    return;
end

while fillingModel.t < 120
    mOld = runTank.mTotalNitrous;
   fillingModel.advanceTime();
   fprintf(['Time: ',num2str(fillingModel.t),'sec, liqHeight: ',num2str(runTank.liquidHeight),', run tank temp: ',num2str(runTank.temp),'K, run tank nitrous: ',num2str(runTank.mTotalNitrous),'Kg, dm: ',num2str(runTank.mTotalNitrous-mOld),'Kg\n']);
    %pause(1);
end

function [P] = atmosphericPressure()
    P = 101325;
end

function [Po,To,fluidType] = intTankConditions()
    global fillSimRunTank;
    global holeInInternalTankHeight;
    if fillSimRunTank.mAir > 0
        fluidType = FluidType.AIR;
    else
        if fillSimRunTank.liquidHeight >= holeInInternalTankHeight
           fluidType = FluidType.NITROUS_LIQUID;
        else
           fluidType = FluidType.NITROUS_GAS;
        end
    end
    
    Po = fillSimRunTank.vapourPressure; %Not QUITE correct (Should use hydrostatic to get) but close enough, and I'm being lazy
    To = fillSimRunTank.temp;
end

function [P] = intTankInletPressure()
    global fillSimRunTank;
    P = fillSimRunTank.pressureAtBaseOfTank;
end

function [Po,To,fluidType] = extTankConditions()
    global fillSimExtTank;
    Po = fillSimExtTank.pressureAtBaseOfTank;
    To = fillSimExtTank.temp;
    fluidType = FluidType.NITROUS_LIQUID;
end
