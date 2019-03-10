%Function to generate a size for an oxidiser tank for a given take off
%temperature (K), required liquid oxidiser mass for flight (Kg), given
%tank cross section area (m^2), given tank ullage fraction (I recommend 20%, 0.2. This is the fraction of vapour when full) and given min excess oxidiser safety margin
%(% of total required as excess). If you see warnings about negative vapour
%mass, INCREASE ullage
%Will return tankHeight in metres
function [tankHeight] = sizeMeATank(tempTakeoff,oxidiserRequired,tankCrossSectionA,tankUllage,minSafetyMargin)
    PVap = SaturatedNitrous.getVapourPressure(tempTakeoff); %Vapour pressure at this temp
    liqDensity = FluidType.NITROUS_LIQUID.getDensity(tempTakeoff,PVap);
    VReq = (oxidiserRequired / liqDensity);
    hReq1 = VReq / tankCrossSectionA;
    tankHeight = (hReq1 / (1-tankUllage)) * 1.15;
    errFunc = @(X) getExcess(X) - (minSafetyMargin*oxidiserRequired);
    tankHeight = abs(real(betterfzero(errFunc,tankHeight,tankHeight*0.9,tankHeight*1.2)));

    function excess = getExcess(tankHeight)
        tankHeight = abs(real(tankHeight));
        tank = tankGen(tempTakeoff,tankUllage,tankHeight,tankCrossSectionA);
        if tank.mVapour < 0
            disp(tank);
           error('Insufficient ullage!?!??'); 
        end
        tank.drainAmountOfLiquid(oxidiserRequired); %Computationally expensive, don't call too often
        excess = tank.mLiquid; %Will go negative once runs out
        margin = excess / oxidiserRequired;
        disp("Latest eval: h: "+tankHeight+", margin: "+margin);
        drawnow;
    end
end