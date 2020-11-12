%Estimates some "real gas" behaviours based on a method suggested in
%"Real_Gas_Thermodynamics", was used before CoolProp/NIST were discovered to exist. 
%Just use CoolProp for values instead of this...


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
            
%             pDelta = PinitGuess-Po;
%             if ~exist('incremSize','var')
%                incremSize = min(100000,abs(pDelta / 10)); %Min of 10 steps, default step size of 100kPa
%                if pDelta < 0 %If need to be decreasing pressure
%                    incremSize = -incremSize; %Should increment by negative amount
%                end
%             end
            
%             PCurrent = Po;
%             TCurrent = To;
%             MCurrent = 0;
%             
%             %Decrement PCurrent until get to faster than required mach
%             %number
%             while MCurrent < MTarget
%                 PGuess = PCurrent + incremSize; %Move P closer to where we want to calculate for
%                 [TCalc,vCalc,M] = RealFlow.getIsentropicTempVelocity(PCurrent,TCurrent,To,CpStagnation,PGuess,fluidType);
%                 
%                 %Update known point that calculations are done relative to
%                 PCurrent = PGuess;
%                 TCurrent = TCalc;
%                 MCurrent = M;
%             end 
            PCurrent = PinitGuess;

            %Our P is currently lower than the correct one, solve for
            %POSITIVE pressure increase that gives correct mach number
            PReq = PCurrent + abs(real(fzero(@errInGuess,0)));
            
            [TCalc,vCalc,MActual] = RealFlow.getIsentropicTempVelocity(Po,To,To,CpStagnation,PReq,fluidType);
            %disp("P FOR M ACTUAL: "+PReq);
            
            %Function to calculate the err in mach number for a
            %desired P given
            function err = errInGuess(PGuessingOffset)
                PGuessingOffset = abs(real(PGuessingOffset)); %Positive real only
                PGuessing = PCurrent + PGuessingOffset;
                [TCalc,vCalc,M] = RealFlow.getIsentropicTempVelocity(Po,To,To,CpStagnation,PGuessing,fluidType);
                err = M - MTarget;
            end
        end
%         
%         %Function to get the temperature of a flow from a known pressure
%         %and temperature that is now at a new known pressure. Should be a
%         %fairly small difference since this assumes the isentropic exponent
%         %is locally constant - which is a good assumption as long as the
%         %pressure gradient you ask about isn't huge.
%         function [TCalc] = getIsentropicTempSmallIncrem(knownPressure,knownTemperature,P,fluidType)
%             PRatio = P / knownPressure;
%             %Isentropic exponent needs a temp and pressure to calculate,
%             %but changes very little with temp and pressure
%             isentropicExponent = fluidType.getRealIsentropicTempPressureExponent(knownTemperature,knownPressure);
%             TCalc = knownTemperature * PRatio^((isentropicExponent-1) / isentropicExponent);
%         end
        
%         %Function to get the temperature and velocity of a flow from a known pressure
%         %and temperature that is now at a new known pressure. Should be a
%         %fairly small difference since this assumes the isentropic exponent
%         %is locally constant - which is a good assumption as long as the
%         %pressure gradient you ask about isn't huge.
%         function [TCalc,vCalc] = getIsentropicTempVelocitySmallIncrem(knownPressure,knownTemperature,To,CpStagnation,P,fluidType)
%             TCalc = RealFlow.getIsentropicTempSmallIncrem(knownPressure,knownTemperature,P,fluidType);
%             vCalc = RealFlow.getIsentropicVelocity(To,TCalc,P,fluidType);
%         end
        
%         %Function to get the velocity and mach number of the isentropic flow for a
%         %given temperature relative to a given total temperature
%         function [vCalc,M] = getIsentropicVelocity(To,T,P,fluidType)
%             a = fluidType.getSpeedOfSound(T,P);
%             %Enthalpy, as implemented now, is only a function of
%             %temperature so give invalid value for pressure.
%             h0 = fluidType.getEnthalpy(To,-1);
%             h = fluidType.getEnthalpy(T,-1);
%             vCalc = sqrt(2*(h0-h));
%             M = vCalc / a;
%             
% %             GPv = fluidType.getRealIsentropicPressureVolumeExponent(T,P);
% %             GTv = fluidType.getRealIsentropicTempVolumeExponent(T,P);
% %             
% %             %From paper on real gas thermodynamics
% %             M = sqrt ( (2/(GPv -1)) * ( (To / T)^( (GPv -1 )/(GTv - 1) ) - 1 ) );
% % 
% %             vCalc = real(M*a);
%         end
        %Function to get the velocity and mach number of the adiabatic flow for a
        %given state relative to upstream
        function [vCalc,M] = getAdiabaticVelocity(TUpstream,PUpstream,vUpstream,T,P,fluidType)
            a = fluidType.getSpeedOfSound(T,P);
            h0 = fluidType.getEnthalpy(TUpstream,PUpstream) + 0.5*vUpstream.^2;
            h = fluidType.getEnthalpy(T,P);
            vCalc = sqrt(2*(h0-h));
            M = vCalc / a;
        end
        
        %Function to get the temperature (K) and velocity (m/s) and Mach number at a point in the
        %flow for a given pressure (Pascal), using given point that has known temperature (K), pressure (Pa)
        %WILL ONLY work for a single phase flow. Be careful because a
        %vapour near the saturation line may need modelling via a biphasic
        %model (See SaturatedNitrous)
        function [TCalc,vCalc,M] = getIsentropicTempVelocity(knownPressure,knownTemperature,knownVelocity,P,fluidType,adiabaticEff)
            if ~exist('adiabaticEff','var')
                adiabaticEff = 1;
            end
            %Uses real-gas-thermodynamics instead of coolprop directly to
            %find isentropic temperature downstream. This is because it
            %agrees nicely with coolprop for values away from saturation
            %line, but as approach saturation line the cool prop values
            %become a little questionable
            
            %PRatio = P / knownPressure;
            %Isentropic exponent needs a temp and pressure to calculate,
            %but changes very little with temp and pressure
            %isentropicExponent = fluidType.getRealIsentropicTempPressureExponent(knownTemperature,knownPressure);
            %disp("PR: "+PRatio+" isentropic exponent: "+isentropicExponent);
            %TDownstreamIsentropic = knownTemperature * PRatio^((isentropicExponent-1) / isentropicExponent);
            entropyUpstream = fluidType.getEntropy(knownTemperature,knownPressure);
            TDownstreamIsentropic = fluidType.getTemperatureFromPressureEntropy(P,entropyUpstream);
            %disp("RGT: "+TCalc+", CoolProp: "+TDownstreamIsentropic);
            hIsen = fluidType.getEnthalpy(TDownstreamIsentropic,P);
            hUpstream = fluidType.getEnthalpy(knownTemperature,knownPressure);
            dhIsentropic = hIsen - hUpstream;
            dh = dhIsentropic * adiabaticEff;
            hTrue = hUpstream+dh;
            TDownstream = fluidType.getTemperatureFromPressureEnthalpy(P,hTrue);
            %disp("PVap downstream: "+fluidType.getVapourPressure(TDownstream));
            TCalc = TDownstream;
            [vCalc,M] = RealFlow.getAdiabaticVelocity(knownTemperature,knownPressure,knownVelocity,TCalc,P,fluidType);
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