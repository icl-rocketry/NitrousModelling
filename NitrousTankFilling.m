classdef NitrousTankFilling < handle
    properties %Allowed to be modified
        internalTank;
        externalTank;
        pipeBetweenTanks;
        pipeVentingInternalTank;
        heightOfHoleInInternalTank; %In metres
        ambientTemp; %In K
        heatTransferCoeffInternalTank;
        externalTankTemp; %In K
        internalTankSurfaceArea;
    end
    properties (SetAccess=private) %Not allowed to be externally modified as may result in non-physical results
        t = 0; %Time in seconds
    end
    properties (Dependent) %Dependent properties
        
    end
    properties (Access=private)
        
    end
    properties (Constant)
       g = 9.81;
       dt = 0.1*10^-3; %1 millisecond step time
    end
    
    methods(Access=private)
    end
    
    methods
        function obj = NitrousTankFilling(internalTankInitial,externalTankInitial,pipeBetweenTanks,pipeVentingInternalTank,heightOfHoleInInternalTank,ambientTemp,heatTransferCoeffInternalTank,internalTankSurfaceArea,externalTankTemp) %Constructor for object
            obj.internalTank = internalTankInitial;
            obj.externalTank = externalTankInitial;
            obj.pipeBetweenTanks = pipeBetweenTanks;
            obj.pipeVentingInternalTank = pipeVentingInternalTank;
            obj.heightOfHoleInInternalTank = heightOfHoleInInternalTank;
            obj.ambientTemp = ambientTemp;
            obj.heatTransferCoeffInternalTank = heatTransferCoeffInternalTank;
            obj.internalTankSurfaceArea = internalTankSurfaceArea;
            obj.externalTankTemp = externalTankTemp;
        end
        
        function advanceTime(obj)
            fprintf('Advancing time... \n');
            %TODO: Change pipe cross sections with time
            
            %Force keeping the external tank at the correct temperature,
            %ignoring the energy required for this
            obj.externalTank.forceSetTemp(obj.externalTankTemp);
            
            %Calculate flow from the external tank to the internal tank
            mdotBetweenTanks = obj.pipeBetweenTanks.calcMassFlow();
            tempIntoIntTank = -1;
            velIntoIntTank = -1;
            %[mdotBetweenTanks,tempIntoIntTank,velIntoIntTank] = obj.pipeBetweenTanks.calcOutlet();
            fprintf(['mdot between tanks: ',num2str(mdotBetweenTanks),' This tick: ',num2str(obj.dt*mdotBetweenTanks),'\n']);
            %pause(2);
            obj.externalTank.drainAmountOfLiquid(obj.dt*mdotBetweenTanks); %Remove from external tank the liquid nitrous
            
            %Add this nitrous to the internal tank
            fprintf(['Adding nitrous mass to internal tank with temperature ',num2str(obj.internalTank.temp),', Cp: ',num2str(FluidType.NITROUS_LIQUID.getCp(obj.externalTank.temp,obj.externalTank.pressureAtBaseOfTank)),', flow temp: ',num2str(tempIntoIntTank),', flow velocity: ',num2str(velIntoIntTank),'\n']);
            CpExt = FluidType.NITROUS_LIQUID.getCp(obj.externalTank.temp,obj.externalTank.pressureAtBaseOfTank);
            TExt = obj.externalTankTemp;
            ho = CpExt * TExt;
            obj.internalTank.addAmountOfNitrous(obj.dt*mdotBetweenTanks,obj.dt*mdotBetweenTanks*ho);
            fprintf(['Run tank temperature after transfer: ',num2str(obj.internalTank.temp),'\n']);
            
            
            %Calculate flow from the internal tank venting out of the
            %rocket
            mDotOut = obj.pipeVentingInternalTank.calcMassFlow();
            fprintf(['mdotOut ',num2str(mDotOut),' Kg/sec\n']);
            mOut = obj.dt * mDotOut;
            
            fprintf(['Should drain ',num2str(mOut),' Kg out\n']);
            
            mDrained = 0;
            while mDrained < mOut && obj.internalTank.mTotalTank > 0
                %If contains air, drain air first
               if obj.internalTank.mAir > 0
                  airDrain = min(mOut,obj.internalTank.mAir);
                  obj.internalTank.drainAmountOfAir(airDrain);
                  mDrained = mDrained + airDrain;
               else %Drain some nitrous
                   lh = obj.internalTank.liquidHeight; %Current height of liquid nitrous in tank
                   %If liquid nitrous is above hole height
                   if obj.heightOfHoleInInternalTank <= lh
                       %Drain liquid nitrous
                       maxLiqToDrain = (lh - obj.heightOfHoleInInternalTank) * obj.internalTank.crossSectionArea * obj.internalTank.liquidDensity;
                       liqDrain = min(mOut,maxLiqToDrain);
                       obj.internalTank.drainAmountOfLiquid(liqDrain);
                       mDrained = mDrained + liqDrain;
                   else
                       %Drain nitrous gas
                       vapourDrain = min(mOut,obj.internalTank.mVapour);
                       obj.internalTank.drainAmountOfVapour(vapourDrain);
                       mDrained = mDrained + vapourDrain;
                   end
               end
            end
            
            %fprintf(['New tank temperature: ',num2str(obj.internalTank.temp),'K\n']);
            
            %Calculate heat transferred to the internal tank in this time
            QInternal = obj.dt * obj.internalTankSurfaceArea * obj.heatTransferCoeffInternalTank * (obj.ambientTemp - obj.internalTank.temp);
            fprintf(['Heat transfer with surroundings: ',num2str(QInternal),'\n']);
            obj.internalTank.addHeat(QInternal); %Transfer this heat
            
            %QExternal = obj.dt * obj.externalTankSurfaceArea * obj.heatTransferCoeffInternalTank * (obj.ambientTemp - obj.externalTank.temp);
            
            fprintf(['Run tank temperature: ',num2str(obj.internalTank.temp),'K\n']);
            fprintf(['Ext tank temperature: ',num2str(obj.externalTank.temp),'\n']);
            fprintf(['Run tank E: ',num2str(obj.internalTank.getInternalEnergy()),' J\n']);
            fprintf(['Ext tank E: ',num2str(obj.externalTank.getInternalEnergy()),' J\n']);
            obj.t = obj.t + obj.dt;
        end
    end
    
end
    