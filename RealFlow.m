%Class to get properties of real flow at various points with isentropic
%and adiabatic equations. 
%
%
%By Eddie Brown
classdef RealFlow
    properties(Constant)
        MIN_ALLOWED_TEMP = 100;
        MAX_ALLOWED_TEMP = 380;
        MIN_ALLOWED_VELOCITY = 0;
        MAX_ALLOWED_VELOCITY = 3*10^8;
    end
    methods (Static)       
        %Calculate the estimated pressure required for a given mach number
        %for a given total pressure and temperature. Uses formula from real
        %gas dynamics to predict, using the pressure-volume polytropic
        %exponent for the stagnation conditions (Which theoretically is
        %able to be considered locally constant), the pressure required to
        %achieve the given mach number. This differs from
        %getRequiredPressureForMach(...) where the pressure-volume
        %polytropic exponent is found for the desired mach number condition
        %iteratively (and expensively)
        %Approximately within 5-10% of the more accurate (but significantly
        %more computationally expensive) values given by getIsentropicConditionsForMach
        function [P] = getEstimateOfRequiredPressureForMach(Po,To,MTarget,fluidType)            
            %Use formula from real gas dynamics to generate initial guess,
            %with wrong temp and P
            GPv = fluidType.getRealIsentropicPressureVolumeExponent(To,Po);
            initGuess = Po / (1 + 0.5*(GPv-1)*MTarget^2)^( (GPv / (GPv - 1)) );
            P = initGuess;
        end
        
        %Function to get the temp, velocity and pressure of flow at a required target mach number for a prescribed total pressure (Pascal) and total temperature (K). incremSize is a small
        %pressure (Pascal) relative to your problem. This function works by stepping from the known pressure condition to JUST BEYOND your given mach condition by this amount. It then will apply a numerical solver (fzero) to find close to exactly the condition you want. If incremSize left blank a sensible default will be used 
        function [TCalc,vCalc,PReq,MActual] = getIsentropicConditionsForMach(Po,To,MTarget,fluidType,incremSize)
            CpStagnation = fluidType.getCp(To,Po);
            
            %Use formula from real gas dynamics to generate initial guess,
            %with wrong temp and P
            GPv = fluidType.getRealIsentropicPressureVolumeExponent(To,Po);
            PinitGuess = Po / (1 + 0.5*(GPv-1)*MTarget^2)^( (GPv / (GPv - 1)) ); %Should be a pretty good guess of pressure for choked flow
            
            pDelta = PinitGuess-Po;
            if ~exist('incremSize','var')
               incremSize = min(100000,abs(pDelta / 10)); %Min of 10 steps, default step size of 100kPa
               if pDelta < 0 %If need to be decreasing pressure
                   incremSize = -incremSize; %Should increment by negative amount
               end
            end
            
            PCurrent = Po;
            TCurrent = To;
            MCurrent = 0;
            
            %Decrement PCurrent until get to faster than required mach
            %number
            while MCurrent < MTarget
                PGuess = PCurrent + incremSize; %Move P closer to where we want to calculate for
                
                [TCalc,vCalc] = RealFlow.getIsentropicTempVelocitySmallIncrem(PCurrent,TCurrent,To,CpStagnation,PGuess,fluidType);
                a = NitrousFluid.getGasSpeedOfSound(TCalc,PGuess);
                M = vCalc / a;
                %disp("v: "+vCalc+", T: "+TCalc+" (P="+PGuess+")");
                
                %Update known point that calculations are done relative to
                PCurrent = PGuess;
                TCurrent = TCalc;
                MCurrent = M;
            end 
            
            %Our P is currently lower than the correct one, solve for
            %POSITIVE pressure increase that gives correct mach number
            PReq = PCurrent + abs(real(fzero(@errInGuess,0)));
            [TCalc,vCalc] = RealFlow.getIsentropicTempVelocitySmallIncrem(PCurrent,TCurrent,To,CpStagnation,PReq,fluidType);
            a = NitrousFluid.getGasSpeedOfSound(TCalc,PReq);
            MActual = vCalc / a;
            %disp("P FOR M ACTUAL: "+PReq);
            
            %Function to calculate the err in mach number for a
            %desired P given
            function err = errInGuess(PGuessingOffset)
                PGuessingOffset = abs(real(PGuessingOffset)); %Positive real only
                PGuessing = PCurrent + PGuessingOffset;
                [TCalc,vCalc] = RealFlow.getIsentropicTempVelocitySmallIncrem(PCurrent,TCurrent,To,CpStagnation,PGuessing,fluidType);
                a = NitrousFluid.getGasSpeedOfSound(TCalc,PGuessing);
                M = vCalc / a;
                err = M - MTarget;
            end
        end
        
        %Function to get the temperature of a flow from a known pressure
        %and temperature that is now at a new known pressure. Should be a
        %fairly small difference since this assumes the isentropic exponent
        %is locally constant - which is a good assumption as long as the
        %pressure gradient you ask about isn't huge.
        function [TCalc] = getIsentropicTempSmallIncrem(knownPressure,knownTemperature,P,fluidType)
            PRatio = P / knownPressure;
            %Isentropic exponent needs a temp and pressure to calculate,
            %but changes very little with temp and pressure
            isentropicExponent = fluidType.getRealIsentropicTempPressureExponent(knownTemperature,knownPressure);
            TCalc = knownTemperature * PRatio^((isentropicExponent-1) / isentropicExponent);
        end
        
        %Function to get the temperature and velocity of a flow from a known pressure
        %and temperature that is now at a new known pressure. Should be a
        %fairly small difference since this assumes the isentropic exponent
        %is locally constant - which is a good assumption as long as the
        %pressure gradient you ask about isn't huge.
        function [TCalc,vCalc] = getIsentropicTempVelocitySmallIncrem(knownPressure,knownTemperature,To,CpStagnation,P,fluidType)
            TCalc = RealFlow.getIsentropicTempSmallIncrem(knownPressure,knownTemperature,P,fluidType);
            vCalc = RealFlow.getIsentropicVelocity(To,TCalc,P,fluidType);
        end
        
        %Function to get the velocity and mach number of the isentropic flow for a
        %given temperature relative to a given total temperature
        function [vCalc,M] = getIsentropicVelocity(To,T,P,fluidType)
            a = fluidType.getSpeedOfSound(T,P);
            GPv = fluidType.getRealIsentropicPressureVolumeExponent(T,P);
            GTv = fluidType.getRealIsentropicTempVolumeExponent(T,P);
            
            %From paper on real gas thermodynamics
            M = sqrt ( (2/(GPv -1)) * ( (To / T)^( (GPv -1 )/(GTv - 1) ) - 1 ) );

            vCalc = real(M*a);
        end
        
        %Function to get the temperature (K) and velocity (m/s) and Mach number at a point in the
        %flow for a given pressure (Pascal), using given point that has known temperature (K), pressure (Pa)
        %and the SAME total temperature To (K). incremSize is a small
        %pressure (Pascal) relative to your problem. This function works by stepping from the known pressure condition to your given condition by this amount. If it left blank a default will be used 
        function [TCalc,vCalc,M] = getIsentropicTempVelocity(knownPressure,knownTemperature,To,CpStagnation,P,fluidType,incremSize)
            if To < knownTemperature
               error('T known > To known??'); 
            end
            pDelta = P-knownPressure;
            if ~exist('incremSize','var')
               incremSize = min(100000,abs(pDelta / 10)); %Min of 10 steps, default step size of 100kPa
               if pDelta < 0 %If need to be decreasing pressure
                   incremSize = -incremSize; %Should increment by negative amount
               end
            end
            
            PCurrent = knownPressure;
            TCurrent = knownTemperature;
            
            while true
                PGuess = PCurrent + incremSize; %Move P closer to where we want to calculate for
                if abs(PGuess - P) <= abs(incremSize) %If this step brought us to or past our pressure we want
                    PGuess = P;
                end
                
                [TCalc,vCalc] = RealFlow.getIsentropicTempVelocitySmallIncrem(PCurrent,TCurrent,To,CpStagnation,PGuess,fluidType);
                %disp("v: "+vCalc+", T: "+TCalc+" (P="+PGuess+")");
                
                if PGuess == P
                    a = NitrousFluid.getGasSpeedOfSound(TCalc,P);
                    M = vCalc / a;
                    drawnow;
                   return; %Values calculated for the correct pressure 
                end
                
                %Update known point that calculations are done relative to
                PCurrent = PGuess;
                TCurrent = TCalc;
            end 
        end
        
        %Function to get the velocity and temperature of flow downstream of
        %a choke point To is in Kelvin, P in pascal, chokedDensity in
        %Kg/m^3, chokedVelocity in m/s, cross section areas in m^3 and fluid type should be of
        %type FluidType (Eg. FluidType.NITROUS_GAS)
        function [TCalc,vCalc] = getVelocityTempForKnownMassFlow(To,P,chokedDensity,chokedVelocity,chokedCrossSectionA,crossSectionA,CpStagnation,fluidType)
             TGuess = P / (fluidType.getGasConstant() * chokedDensity); %Estimate from eqn of state with compressibility of 1
             solution = real(betterfzero(@errFunc,TGuess,180,380));
             TCalc = solution(1);
             vCalc = getIsentropicVelocity(To,TCalc,P,fluidType);
             
             %fprintf(['T ',num2str(TCalc),' P ',num2str(P),' M ','choked',' vel ',num2str(vCalc),'\n']);
             
             %Function for solver to use to find correct T
             function err = errFunc(T)
                 T = real(T);
                 %From conservation of mass, and To equation
                 vel1 = (RealFlow.getIsentropicVelocity(To,T,P,fluidType));
                 rho = fluidType.getDensity(T,P);
                 %From conservation of mass with changing cross section area
                 %and density
                 vel2 = (chokedDensity / rho) * chokedVelocity * (chokedCrossSectionA/crossSectionA);
                 err = vel2 - vel1;
                 %fprintf(['T: ',num2str(T),' vel err ',num2str(err),' vel ',num2str(vel1),'\n']);
             end
        end
    end
    
end