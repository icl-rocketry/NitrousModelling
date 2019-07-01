%Enumeration to represent types of fluids
classdef FluidType
   enumeration
      AIR(@(T,P) 1005,@(T,P) 718, @(T,P) P/(287*T),287, @(T,P) sqrt(1.401*287*T), @(P) error('Unsupported'), @(T,P) 1005 * T, @(T,P) error('Unsupported'),@(T,P) error('Unsupported'),@(T,s) error('Unsupported'),@(P,s) error('Unsupported'),@(T,s) error('Unsupported'),@(P,s) error('Unsupported')), 
      NITROUS_LIQUID(@NitrousFluid.getLiquidCp,@NitrousFluid.getLiquidCv,@NitrousFluid.getLiquidDensity,@NitrousFluid.getGasConstant,@NitrousFluid.getLiquidSpeedOfSound,@SaturatedNitrous.getVapourPressure,@NitrousFluid.getLiquidSpecificEnthalpy, @NitrousFluid.getLiquidIsobaricCoeffOfExpansion,@NitrousFluid.getSaturatedLiquidEntropy,@NitrousFluid.getSaturatedLiquidPressureFromEntropy,@NitrousFluid.getSaturatedLiquidTemperatureFromEntropy,@NitrousFluid.getLiquidPressureForTemperatureEnthalpy,@NitrousFluid.getLiquidTemperatureForPressureEnthalpy), 
      NITROUS_GAS(@NitrousFluid.getGasCp,@NitrousFluid.getGasCv,@NitrousFluid.getGasDensity,@NitrousFluid.getGasConstant,@NitrousFluid.getGasSpeedOfSound,@SaturatedNitrous.getVapourPressure,@NitrousFluid.getGasSpecificEnthalpy, @NitrousFluid.getGasIsobaricCoeffOfExpansion,@NitrousFluid.getGasEntropy,@NitrousFluid.getGasPressureForTemperatureEntropy,@NitrousFluid.getGasTemperatureForPressureEntropy,@NitrousFluid.getGasPressureForTemperatureEnthalpy,@NitrousFluid.getGasTemperatureForPressureEnthalpy);
   end
   properties
      CpHandle;
      CvHandle;
      enthalpyHandle;
      entropyHandle;
      pressureFromTempEntropyHandle;
      temperatureFromPressureEntropyHandle;
      pressureFromTempEnthalpyHandle;
      temperatureFromPressureEnthalpyHandle;
      densityHandle;
      gasConstant;
      speedOfSoundHandle;
      vapourPressureHandle;
      isobaricCoeffOfExpansionHandle;
   end
   methods
       function obj = FluidType(CpHandle,CvHandle,densityHandle,gasConstant,speedOfSoundHandle,vapourPressureHandle,enthalpyHandle,isobaricCoeffOfExpansionHandle,entropyHandle,pressureFromTempEntropyHandle,temperatureFromPressureEntropyHandle,pressureFromTempEnthalpyHandle,temperatureFromPressureEnthalpyHandle)
         obj.CpHandle = CpHandle;
         obj.CvHandle = CvHandle;
         obj.densityHandle = densityHandle;
         obj.gasConstant = gasConstant;
         obj.speedOfSoundHandle = speedOfSoundHandle;
         obj.vapourPressureHandle = vapourPressureHandle;
         obj.enthalpyHandle = enthalpyHandle;
         obj.isobaricCoeffOfExpansionHandle = isobaricCoeffOfExpansionHandle;
         obj.entropyHandle = entropyHandle;
         obj.pressureFromTempEntropyHandle = pressureFromTempEntropyHandle;
         obj.temperatureFromPressureEntropyHandle = temperatureFromPressureEntropyHandle;
         obj.pressureFromTempEnthalpyHandle = pressureFromTempEnthalpyHandle;
         obj.temperatureFromPressureEnthalpyHandle = temperatureFromPressureEnthalpyHandle;
       end
       
       function P = getVapourPressure(obj,T)
            P = obj.vapourPressureHandle(T);
       end
       
       %Define fluidType.getCp(T,P) so simpler to get Cp of a given fluid,
       %don't have to check each fluid type and call correct method as it's
       %defined within the enum
       function Cp = getCp(obj,T,P)
           if length(T) > 1 && length(P) == 1
               for i=1:length(T)
                   Cp(i) = obj.CpHandle(T(i),P);
               end
           elseif length(T) == 1 && length(P) > 1
               for i=1:length(P)
                   Cp(i) = obj.CpHandle(T,P(i));
               end
           elseif length(T) == length(P)
               for i=1:length(T)
                   Cp(i) = obj.CpHandle(T(i),P);
               end
           end
       end
       
       %Get the specific internal energy of the fluid using u = h-pv and
       %enthalpy data for the given fluid
       function u = getSpecificInternalEnergy(obj,T,P)
            h = obj.getSpecificEnthalpy(T,P);
            v = 1 / obj.getDensity(T,P);
            %H = U + PV
            u = h - P*v;
       end
       
       function h = getSpecificEnthalpy(obj,T,P)
            h = obj.enthalpyHandle(T,P);
%            v = 1 / obj.getDensity(T,P);
%            h = obj.getCv(T,P) * T + P * v;
       end
       
       %Function to predict enthalpy at a given temp and pressure of the
       %fluid given a known enthalpy at a known temp and pressure. 
       %Uses dh = CpdT + [v(1-aT)]dP
       function val = predictEnthalpy(obj,knownH,knownT,knownP,T,P)
            dP = P-knownP; %Difference in pressure
            stepSize = 10000; %Step size to use for numerical integration
            steps = ceil(abs(dP / stepSize)); %Number of discrete steps for integration
            
            integrand = @(Pi) (1/obj.getDensity(knownT,Pi))*(1-(obj.isobaricCoeffOfExpansionHandle(T,Pi))*T);
            dh1 = legendreIntegral(integrand,max(10,steps),SaturatedNitrous.getVapourPressure(T),P);
            
            %Numerically integrate CpdT at constant P
            dT = T-knownT;
            stepSize = 10000;
            steps = ceil(abs(dT / stepSize));
            
            integrand2 = @(Ti) obj.getCp(Ti,knownP);
            dh2 = legendreIntegral(integrand2,max(10,steps),knownT,T);
            %fprintf(['dh: ',num2str(dh),'\n']);
            
            val = knownH + dh1 + dh2; 
       end
       
       %Get the Cv of the gas using method supplied by gas
       function Cv = getCv(obj,T,P)
            Cv = obj.CvHandle(T,P);
       end
       
       %Function to calculate real Cv of ANY given gas from real Cp and density
       %temperature. Commented out as Cv handle specifies for each gas (So
       %more accurate data can be supplied if given gas has more accurate
       %data)
%        function Cv = calcRealCv(obj,T,P)
%            smallIncrem = 1*10^-6;
%            %Density and V at this point
%            V0 = 1/obj.getDensity(T-smallIncrem,P);
%            V1 = 1/obj.getDensity(T,P-smallIncrem);
%            
%            %Density with small increm in temp and with small increm in P
%            %respectively
%            V2 = 1/obj.getDensity(T+smallIncrem,P);
%            V3 = 1/obj.getDensity(T,P+smallIncrem);
%            
%            %Central finite difference approximation           
%            %Partial dV/dT at constant pressure
%            dVdT = (V2 - V0) / (2*smallIncrem);
%            %Partial dV/dP at constant temperature
%            dVdP = (V3 - V1) / (2*smallIncrem);
%            
%            Cv = obj.getCp(T,P) + (T) * ( (dVdT)^2 / (dVdP) );
%        end
       
       function gamma = getGamma(obj,T,P)
           gamma = obj.getCp(T,P) / obj.getCv(T,P);
       end
       
       %Gets speed of sound of fluid in m/s using method determined for
       %this fluid (Eg. for nitrous 1/(rho*betaS))
       function a = getSpeedOfSound(obj,T,P)
          a = obj.speedOfSoundHandle(T,P); 
       end
       
       %Gets speed of sound of fluid in m/s using real gas isentropic
       %pressure volume exponent. Should use getSpeedOfSound for more
       %accurate results as data derived
       function a = getSpeedOfSound2(obj,T,P)
           Z = obj.getCompressibilityFactor(T,P);
           R = obj.getGasConstant();
          a = sqrt(obj.getRealIsentropicPressureVolumeExponent(T,P) * Z * R * T); 
       end
%        
%        %Get speed of sound of fluid based on real gas thermodynamics
%        %bulk modulus approach
%        function c = getSpeedOfSound3(obj,T,P)
%            K = obj.getBulkModulus2(T,P);
%            density = obj.getDensity(T,P);
%            c = sqrt(K / density);
%        end
       
       %Get the isothermal compressibility of the fluid based on real gas
       %thermodynamics paper approach
       function comp = getIsothermalCompressibility2(obj,T,P)
           gammaPT = obj.getRealIsentropicTempPressureExponent(T,P);
           Cp = obj.getCp(T,P);
           v = 1 / obj.getDensity(T,P);
           comp = ((gammaPT-1)/gammaPT) * (Cp / (P*v));
       end
       
       %Get bulk modulus of fluid based on real gas thermodynamics approach
       function K = getBulkModulus2(obj,T,P)
           gammaPV = obj.getRealIsentropicPressureVolumeExponent(T,P);
           K = gammaPV * P;
       end
       
       function Z = getCompressibilityFactor(obj,T,P)
           rho = obj.getDensity(T,P);
           R = obj.getGasConstant();
            Z = P / (rho * R * T); 
       end
       
       function rho = getDensity(obj,T,P)
           if length(T) > 1 && length(P) == 1
               for i=1:length(T)
                   rho(i) = obj.densityHandle(T(i),P);
               end
           elseif length(T) == 1 && length(P) > 1
               for i=1:length(P)
                   rho(i) = obj.densityHandle(T,P(i));
               end
           elseif length(T) == length(P)
               for i=1:length(T)
                   rho(i) = obj.densityHandle(T(i),P);
               end
           end
       end
       
       function R = getGasConstant(obj)
           R = obj.gasConstant();
       end
       
       %Bulk modulus in Pascal
       function val = getBulkModulus(obj,T,P)
           val = obj.getDensity(T,P) * (obj.getSpeedOfSound(T,P)^2);
       end
       
       % Equivalent to isentropic exponent, just less accurate as higher
       % uncertainty
%        function val = findPolytropicTempPressureExponent(obj,refTemp,refPressure,T,P)
%            smallChange = 1*10^-6;
%            KP2 = obj.getBulkModulus(T,P+0.5*smallChange);
%            KP1 = obj.getBulkModulus(T,P-0.5*smallChange);
%            dKdP = (KP2 - KP1) / smallChange; %partial dK dP
%            
%            KT2 = obj.getBulkModulus(T+0.5*smallChange,P);
%            KT1 = obj.getBulkModulus(T-0.5*smallChange,P);
%            dKdT = (KT2 - KT1) / smallChange; %partial dK dP
%            
%            val = fzero(@errFunc,obj.getRealIsentropicTempPressureExponent(T,P));
%            
%            function err = errFunc(a)
%                A = refTemp / refPressure^((a-1)/a);
%                B = (a-1)/a;
%                C = P^(-1/a);
%                err = dKdP + dKdT * (A*B*C) - a;
%            end
%        end

%         %Uses real gas thermodynamics to calculate the entropy change
%         %from state T1, P1 to state T2, P2. See real gas thermodynamics
%         %paper for derivation. Should be a small difference since Cp and
%         %pressure temperature isentropic exponent assumed locally constant
%         function ds = calcEntropyChangeForSmallDifference(obj,T1,T2,P1,P2)
%             Cp = obj.getCp(0.5*(T1+T2),0.5*(P1+P2));
%             gammaPT = obj.getRealIsentropicTempPressureExponent(0.5*(T1+T2),0.5*(P1+P2));
%             ds = Cp*log(T2/T1) - Cp*((gammaPT-1)/gammaPT)*log(P2/P1);
%         end
        
        %Uses real gas thermodynamics equation to calculate the entropy
        %change from state T1, P1 to state T2, P2
        function ds = calcEntropyChange(obj,T1,T2,P1,P2)
            %TODO Re-write term 1 and term 2 as functions that calculate
            %for each T and each P as params passed in are arrays
            for i=1:length(T1)
                term1 = @(T) obj.getCp(T,P1(i)) .* (1./T);
                GPT = @(P) obj.getRealIsentropicTempPressureExponent(T1(i),P);
                term2 = @(P) obj.getCp(T1(i),P) .* ((GPT(P)-1)./GPT(P)) .* (1./P);
                
                ds(i) = integral(term1,T1(i),T2(i)) - integral(term2,P1(i),P2(i));
            end
        end
        
        function s = getEntropy(obj,T,P)
            try 
                s = obj.entropyHandle(T,P);
            catch
                s = obj.entropyHandle(T);
            end
        end
        
        function P = getPressureFromTemperatureEntropy(obj,T,s)
            try 
                P = obj.pressureFromTempEntropyHandle(T,s); 
            catch
                P = obj.pressureFromTempEntropyHandle(s); 
            end
        end
        
        function T = getTemperatureFromPressureEntropy(obj,P,s)
            try 
                T = obj.temperatureFromPressureEntropyHandle(P,s); 
            catch
                T = obj.temperatureFromPressureEntropyHandle(s);  
            end
        end
        
        function P = getPressureFromTemperatureEnthalpy(obj,T,h)
            P = obj.pressureFromTempEnthalpyHandle(T,h);
        end
        
        function P = getTemperatureFromPressureEnthalpy(obj,P,h)
            P = obj.temperatureFromPressureEnthalpyHandle(P,h);
        end
        
%         %DO NOT USE. uses calcEntropyChange and NIST data as a
%         %reference point to determine absolute entropy values for modelling
%         %saturated pipe flow
%         function s = getEntropyTest(obj,T,P)
%            %Entropy of nitrous liquid relative to reference
% %            refT = 182.33;
% %            refP = obj.getVapourPressure(182.33);
%            refT = 150;
%            refP = 101.3*10^3;
%            %s = obj.calcEntropyChange(refT,T,refP,P) + 2550; %2550 is a constant for nitrous liquid
%            s = obj.calcEntropyChange(refT,T,refP,P) + 4445.454545; %2550 is a constant for nitrous liquid
%         end
       
       %See 'real gas thermodynamics' paper for derivation
        function val = getRealIsentropicTempPressureExponent(obj,T,P)
            R = obj.getGasConstant();
            Cp = obj.getCp(T,P);
            smallIncrem = 1*10^-7;
            Z1 = obj.getCompressibilityFactor(T-smallIncrem,P);
            Z = obj.getCompressibilityFactor(T,P);
            Z2 = obj.getCompressibilityFactor(T+smallIncrem,P);
            dZdT = (Z2-Z1) ./ (2*smallIncrem); %Partial derivative of dZ/dT with P constant
            val = (1-(R./Cp).*(Z + T.*dZdT)  ).^-1;
        end
        
        %See 'real gas thermodynamics' paper for derivation
        function val = getRealIsentropicTempVolumeExponent(obj,T,P)
            R = obj.getGasConstant();
            Cv = obj.getCv(T,P);
            Z = obj.getCompressibilityFactor(T,P);
            rho = P / (Z * R * T); %Density to keep constant
            smallIncrem = 1*10^-7; %Small increm in T for calculating dZ/dT with density constant
            T2 = T+smallIncrem;
            T0 = T-smallIncrem;
            %Find P that keeps rho the same for this change in T, with P/2
            %and 3P/2 as min and max bounds
            P2 = betterfzero(@errDensity,P,P*0.5,P*1.5);
            P2 = real(P2);
            Z2 = obj.getCompressibilityFactor(T2,P2);
            P0 = betterfzero(@errDensity2,P,P*0.5,P*1.5);
            P0 = real(P0);
            Z0 = obj.getCompressibilityFactor(T0,P0);
            dZdT = (Z2-Z0) / (2*smallIncrem); %Partial derivative of dZ/dT with vol/density constant
            
            val = 1 + (R/Cv) * (Z + T*dZdT); %Formula from paper
            
            function err = errDensity(P) %Error in density
                P = real(P);
                Z2 = obj.getCompressibilityFactor(T2,P);
                rho2 = P / (Z2 * R * T2);
                err = rho2 - rho;
            end
            
            function err = errDensity2(P) %Error in density
                P = real(P);
                Z0 = obj.getCompressibilityFactor(T0,P);
                rho2 = P / (Z0 * R * T0);
                err = rho2 - rho;
            end
        end
        
        %See 'real gas thermodynamics' paper for derivation
        function val = getRealIsentropicPressureVolumeExponent(obj,T,P)
            Tv = obj.getRealIsentropicTempVolumeExponent(T,P);
            PT = obj.getRealIsentropicTempPressureExponent(T,P);
            val = (Tv -1) * (PT / (PT-1));
        end
   end
end