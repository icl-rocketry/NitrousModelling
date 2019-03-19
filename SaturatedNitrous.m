%Class to get properties of nitrous in saturated state (Liquid and gas
%inside a contained above boiling temperature)
%Valid UP TO the critical temperature only (309.57K)
%Uses equations from http://edge.rit.edu/edge/P07106/public/Nox.pdf
%for modelling properties
%Very similar can be seem implemented in C++ in the following paper:
%http://www.aspirespace.org.uk/downloads/Modelling%20the%20nitrous%20run%20tank%20emptying.pdf
%By Eddie Brown
classdef SaturatedNitrous
    properties(Constant)
        T_CRIT = 309.57; %Critical temp of nitrous in kelvin
        P_CRIT = 72.51 * 10^5; %Critical pressure of nitrous, Pa
        DENSITY_CRIT = 452; %Critical density of nitrous, Kg/m^3
    end
    methods (Static) %Static methods: Methods that don't need object instance
        %Get's the specific entropy of saturated nitrous flow mixture of
        %given quality X, temperature T and pressure P
        function s = getSpecificEntropy(X,T,P)
            sVapour = FluidType.NITROUS_GAS.getEntropy(T,P);
            sLiquid = FluidType.NITROUS_LIQUID.getEntropy(T,P);
            s = sVapour*X + sLiquid*(1-X);
        end
        
        %Get's the specific enthalpy of saturated nitrous flow mixture of
        %given quality X, temperature T and pressure P
        function h = getSpecificEnthalpy(X,T,P)
            hVapour = FluidType.NITROUS_GAS.getSpecificEnthalpy(T,P);
            hLiquid = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(T,P);
            h = hVapour*X + hLiquid*(1-X);
        end
        
        %Bubble growth time characteristic as defined in "Modeling feed
        %system flow physics for self pressurising propellants" by Jonny
        %Dyer https://arc.aiaa.org/doi/pdf/10.2514/6.2007-5702. 1 denotes
        %upstream, 2 denotes downstream. Units are SI
        function Tb = getBubbleGrowthTimeCharacteristic(P1,T1,P2)
            rhoLiq = FluidType.NITROUS_LIQUID.getDensity(T1,P1);
            vapourPressure = SaturatedNitrous.getVapourPressure(T1);
            Tb = sqrt(1.5 * (rhoLiq/(vapourPressure-P2)));
        end
        
        %Residence time (proportional to how long fluid will take to do this) characteristic as defined in "Modeling feed
        %system flow physics for self pressurising propellants" by Jonny
        %Dyer https://arc.aiaa.org/doi/pdf/10.2514/6.2007-5702. 1 denotes
        %upstream, 2 denotes downstream. Units are SI. Length is distance
        %in metres flow is travelling over for this pressure change
        function Tr = getResidenceTime(P1,T1,P2,length)
            rhoLiq = FluidType.NITROUS_LIQUID.getDensity(T1,P1);
            Tr = length*sqrt(rhoLiq/(2*(P1-P2)));
        end
        
        %Get the non-equilibrium flow parameter as defined in "Modeling feed
        %system flow physics for self pressurising propellants" by Jonny
        %Dyer https://arc.aiaa.org/doi/pdf/10.2514/6.2007-5702. 1 denotes
        %upstream, 2 denotes downstream. Units are SI. Length is distance
        %flow is travelling over for this pressure change. If results given
        %by this are not as desired, try adjusting "length" until it
        %changes the output by the correct amount. (If normal flow
        %coefficient applied to mass flow isn't enough)
        %Non-equilibrium flow parameter is a ratio of bubble time to
        %residence time (If large then bubbles take longer to form relative
        %to time taken for this flow and therefore have little time to form, 
        %if small then bubbles have a long time to form
        function k = getNonEquilibriumFlowParameter(P1,T1,P2,length)
            k = SaturatedNitrous.getBubbleGrowthTimeCharacteristic(P1,T1,P2) / SaturatedNitrous.getResidenceTime(P1,T1,P2,length);
        end
        
        %Get the downstream mass flow rate predicted by an ideal incompressible
        %assumption for a given upstream pressure, temperature and velocity
        %and given downstream pressure. G is mass flow per unit area, all
        %units SI. Temperature given by this also assumes flow is
        %adiabatic and quality given by this assumes isentropic. Set ignoreTempAndQuality to true if you don't want to calculate
        %the temperature and quality downstream (will make this 1000x faster to
        %compute)
        function [X2,T2,v2,h2,G] = getDownstreamIncompressibleMassFlowTemp(X1,P1,T1,P2,v1,ignoreTempAndQuality)
            rho = FluidType.NITROUS_LIQUID.getDensity(T1,P1); %Upstream density, as incompressible assume is constant
            G = sqrt(2.*rho.*(P1-P2) + rho.^2.*v1.^2);
            v2 = G ./ rho; %G = rho * v
            %Assume total enthalpy constant to figure out T
            h0 = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(T1,P1) + (v1.^2)./2;
            h2 = h0 - (v2.^2/2); %Downstream enthalpy
            
            if ~exist('ignoreTemp','var') || ~ignoreTempAndQuality
                T2 = real(fzero(@T2Err,T1));
                s1 = SaturatedNitrous.getSpecificEntropy(X1,T1,P1);
                X2 = real(fzero(@(X) real(SaturatedNitrous.getSpecificEntropy(real(X),T2,P2)-s1),X1)); %Find quality downstream that satisfies isentropic condition
            else
                T2 = -1; 
                X2 = 0;
            end
            function err = T2Err(TGuess) %Error in guess for downstream temperature (to match the enthalpy)
                h =  FluidType.NITROUS_LIQUID.getSpecificEnthalpy(real(TGuess),P2);
                err = h-h2; %0 if downstream enthalpy matches
            end
        end
        
        %Get the flow rate and characteristics predicted by a
        %non-homogenous non equilibrium model as defined By Dyer in
        %"Modeling feed system physics for self pressurising propellants"
        %(https://arc.aiaa.org/doi/pdf/10.2514/6.2007-5702)
        %and with the correction to the model coefficients made as given by
        %Brian Solomon in "Engineering Model to calculate mass flow rate of
        %a two phase saturated fluid through an injector orifice"
        %(https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=1110&context=gradreports)
        %Some slight extra handling added to give good answers even if
        %input not pefectly on saturation to start
        %Mass flow answers given by this should be corrected with an empirically
        %determined coefficient
        function [X2,T2,v2,h2,G] = getDownstreamSaturatedNHNEFlowCond(X1,T1,P1,P2,v1,length)
            %While P1 larger than vapour pressure for given temperature by
            %more than 0.1%, and quality close enough to liquid, drop
            %pressure isentropically until reach saturation line
            while ((P1-SaturatedNitrous.getVapourPressure(T1))/SaturatedNitrous.getVapourPressure(T1)) > 0.001 && X1 < 0.01
                %Drop pressure to vapour pressure and keep a track of
                %velocity,temp etc...
                vapourPressure = max(SaturatedNitrous.getVapourPressure(T1),P2);
                h1 = SaturatedNitrous.getSpecificEnthalpy(X1,T1,P1); %Enthalpy upstream
                [T1,~,~] = RealFlow.getIsentropicTempVelocity(P1,T1,T1,FluidType.NITROUS_LIQUID.getCp(T1,P1),vapourPressure,FluidType.NITROUS_LIQUID);
                h2 = SaturatedNitrous.getSpecificEnthalpy(X1,T1,vapourPressure); %Enthalpy after pressure drop
                KEGain = h1-h2; %specific KE gain is difference in fluid enthalpies upstream and downstream
                v1 = sqrt(2*(0.5*v1.^2 + KEGain)); %Energy balance, assuming no other losses
                P1 = vapourPressure;
                if P1 == P2 %Flow never reached saturation line
                    X2 = X1;
                    T2 = T1;
                    v2 = v1;
                    %h2 = h2;
                    rhoLiq = FluidType.NITROUS_LIQUID.getDensity(T2,P2);
                    G = rhoLiq * v2;
                    return;
                end
            end
            k = SaturatedNitrous.getNonEquilibriumFlowParameter(P1,T1,P2,length);
            incompressibleCoeff = (1-1/(1+k)); %From correction to dyer by solomon
            hemCoeff = (1/(1+k)); %From correction to dyer by solomon
            
            %Calculate flow both incompressibly and with HEM
            [X2Inc,T2Inc,v2Inc,h2Inc,GInc] = SaturatedNitrous.getDownstreamIncompressibleMassFlowTemp(X1,P1,T1,P2,v1);
            [X2Hem,T2Hem,v2Hem,h2Hem,GHem] = SaturatedNitrous.getDownstreamIsentropicSaturatedHEMFlowCond(X1,T1,P1,P2,v1);
            %Weighting of flow properties as in Dyer
            X2 = incompressibleCoeff * X2Inc + hemCoeff * X2Hem;
            T2 = incompressibleCoeff * T2Inc + hemCoeff * T2Hem;
            v2 = incompressibleCoeff * v2Inc + hemCoeff * v2Hem;
            h2 = incompressibleCoeff * h2Inc + hemCoeff * h2Hem;
            G = incompressibleCoeff * GInc + hemCoeff * GHem;
        end
        
        %Get downstream conditions of a flow that is modelled by isentropic
        %saturated homogenous equilibrium model - very much an imperfect
        %model, but it seems like better models for this kind of flow don't
        %really exist and best approach is correcting something like this
        %with empirical data. G is mass flow per unit area (Kg/sec per
        %m^2). X is the flow "quality", eg. fraction that is vapour
        function [X2,T2,v2,h2,G] = getDownstreamIsentropicSaturatedHEMFlowCond(X1,T1,P1,P2,v1)
            %While P1 larger than vapour pressure for given temperature by
            %more than 0.1%, and quality close enough to liquid, drop
            %pressure isentropically until reach saturation line
            while ((P1-SaturatedNitrous.getVapourPressure(T1))/SaturatedNitrous.getVapourPressure(T1)) > 0.001 && X1 < 0.01
                %Drop pressure to vapour pressure and keep a track of
                %velocity,temp etc...
                vapourPressure = max(SaturatedNitrous.getVapourPressure(T1),P2);
                h1 = SaturatedNitrous.getSpecificEnthalpy(X1,T1,P1); %Enthalpy upstream
                [T1,~,~] = RealFlow.getIsentropicTempVelocity(P1,T1,T1,FluidType.NITROUS_LIQUID.getCp(T1,P1),vapourPressure,FluidType.NITROUS_LIQUID);
                h2 = SaturatedNitrous.getSpecificEnthalpy(X1,T1,vapourPressure); %Enthalpy after pressure drop
                KEGain = h1-h2; %specific KE gain is difference in fluid enthalpies upstream and downstream
                v1 = sqrt(2*(0.5*v1.^2 + KEGain)); %Energy balance, assuming no other losses
                P1 = vapourPressure;
                if P1 == P2 %Flow never reached saturation line
                    X2 = X1;
                    T2 = T1;
                    v2 = v1;
                    %h2 = h2;
                    rhoLiq = FluidType.NITROUS_LIQUID.getDensity(T2,P2);
                    G = rhoLiq * v2;
                    return;
                end
            end
            
            %Assuming now actually on saturation line
            s1 = SaturatedNitrous.getSpecificEntropy(X1,T1,P1); %Entropy upstream
            h1 = SaturatedNitrous.getSpecificEnthalpy(X1,T1,P1); %Enthalpy upstream
            T2 = SaturatedNitrous.getSaturationTemperature(P2); %Drops to saturation temp of new pressure
            X2 = real(fzero(@(X) real(SaturatedNitrous.getSpecificEntropy(real(X),T2,P2)-s1),X1)); %Find quality downstream that satisfies isentropic condition
            if X2 < 0 || X2 > 1
               error('INVALID x2 encountered'); 
            end
            h2 = SaturatedNitrous.getSpecificEnthalpy(X2,T2,P2); %Gets the specific enthalpy downstream
            KEGain = h1-h2; %specific KE gain is difference in fluid enthalpies upstream and downstream
            v2 = sqrt(2*(0.5*v1.^2 + KEGain)); %Energy balance, assuming no other losses
            rhoLiq = FluidType.NITROUS_LIQUID.getDensity(T2,P2);
            rhoVap = FluidType.NITROUS_GAS.getDensity(T2,P2);
            rhoMixture = (1-X2)*rhoLiq + X2*rhoVap;
            G = rhoMixture * v2;
        end
        
        %Get saturation temperature (K) for a given pressure (Pascal), if one does not
        %exist might return -1 or might return errornous solution
        function T = getSaturationTemperature(P)
            if P > SaturatedNitrous.P_CRIT
               warning('No saturation temperature for given pressure, fluid is supercritical'); 
               T = -1;
               return;
            end
            
            T = fzero(@(T) (real(SaturatedNitrous.getVapourPressure(real(T)) - P)),290);
            if abs(imag(SaturatedNitrous.getVapourPressure(real(T)))) > 0
                warning('No saturation temperature for given pressure, fluid cannot be liquid'); 
               T = -1; %Does not exist 
            end
        end
        
        %Gets the vapour pressure (Pa) of saturated nitrous at a given
        %temperature (K)
        function P = getVapourPressure(T)
            TRel = T / SaturatedNitrous.T_CRIT; %Relative temp
            b = [-6.71893,1.35966,-1.3779,-4.051]; %Constants in equation
            %Equation from paper
            logPressureRatio = (1/TRel) * ( b(1)*(1-TRel) + b(2)*(1-TRel)^1.5 + b(3)*(1-TRel)^2.5 + b(4)*(1-TRel)^5 );
            P = exp(logPressureRatio) * SaturatedNitrous.P_CRIT; %Calculate vapour pressure
        end
        
        %Gets the liquid density (Kg/m^3) of saturated nitrous at a given
        %temperature (K)
        function rho = getLiquidDensity(T)
            TRel = T / SaturatedNitrous.T_CRIT; %Relative temp
            b = [1.72328,-0.83950,0.51060,-0.10412]; %Constants in equation
            %Equation from paper
            logDensityRatio =  b(1)*(1-TRel)^(1/3) + b(2)*(1-TRel)^(2/3) + b(3)*(1-TRel)^(1) + b(4)*(1-TRel)^(4/3);
            rho = exp(logDensityRatio) * SaturatedNitrous.DENSITY_CRIT; %Calculate vapour pressure
        end
        
        %Gets the vapour density (Kg/m^3) of saturated nitrous at a given
        %temperature (K)
        function rho = getVapourDensity(T)
            TRel = SaturatedNitrous.T_CRIT/T; %Inverse of Relative temp
            b = [-1.00900,-6.28792,7.50332,-7.90463,0.629427]; %Constants in equation
            %Equation from paper
            logDensityRatio =  b(1)*(TRel-1)^(1/3) + b(2)*(TRel-1)^(2/3) + b(3)*(TRel-1)^(1) + b(4)*(TRel-1)^(4/3) + b(5)*(TRel-1)^(5/3);
            rho = exp(logDensityRatio) * SaturatedNitrous.DENSITY_CRIT; %Calculate vapour pressure
        end
        
        %Gets the specific enthalpy (J/Kg) of the saturated liquid at a given
        %temperature (K). Valid up to 35C.  The specific enthalpy values are calculated on the basis that the enthalpy of the ideal gas is
%zero at 298.15 K. HAS a weird offset
        function h = getLiquidSpecificEnthalpy(T)
            b = [-200, 116.042,-917.225,794.779,-589.587];
            Tr = T/SaturatedNitrous.T_CRIT;
            h = b(1) + b(2)*(1-Tr)^(1/3) + b(3)*(1-Tr)^(2/3) + b(4)*(1-Tr) + b(5)*(1-Tr)^(4/3);
            h = 1000*h; %kJ -> J
        end
        
        function h = getAbsoluteLiquidSpecificEnthalpy(T)
            %Enthalpy data is offset by 769576.6574
           h = 769576.6574 + SaturatedNitrous.getLiquidSpecificEnthalpy(T);
        end
        
        %Gets the specific enthalpy (J/Kg) of the saturated vapour at a given
        %temperature (K) Valid up to 36C  The specific enthalpy values are calculated on the basis that the enthalpy of the ideal gas is
%zero at 298.15 K. 
        function h = getVapourSpecificEnthalpy(T)
            b = [-200,440.055,-459.701,434.081,-485.338];
            Tr = T/SaturatedNitrous.T_CRIT; %Relative temp
            h = b(1) + b(2)*(1-Tr)^(1/3) + b(3)*(1-Tr)^(2/3) + b(4)*(1-Tr) + b(5)*(1-Tr)^(4/3);
            h = 1000*h; %kJ -> J
        end
        
        function h = getAbsoluteVapourSpecificEnthalpy(T)
            %Enthalpy data is offset by 769576.6574
           h = 769576.6574 + SaturatedNitrous.getVapourSpecificEnthalpy(T);
        end
        
        %Gets the specific enthalpy (J/Kg) of the ideal nitrous oxide gas at a given
        %temperature (K). NOTE this is a bit weird as it's adjusted to be 0
        %at 25C
        function h = getIdealGasSpecificEnthalpy(T)
            b = [-209.559,61.3277,-52.5969,249.352,-38.4368];
            Tr = T/SaturatedNitrous.T_CRIT; %Relative temp
            h = b(1) + b(2)*(Tr)^(1/2) + b(3)*Tr + b(4)*Tr^(3/2) + b(5)*Tr^2;
            h = 1000*h; %kJ -> J
        end
        
        %Gets the liquid enthalpy (latent heat) of vaporisation (J/Kg) of the saturated liquid at a
        %given temperature (K)
        function H = getLiquidEnthalpyOfVaporisation(T)
            H = SaturatedNitrous.getVapourSpecificEnthalpy(T) - SaturatedNitrous.getLiquidSpecificEnthalpy(T);
        end
        
        %Gets the isobaric heat capacity (Cp, in J/KgK) of the saturated liquid at a
        %given temperature (J). Valid up to 30C
        function Cp = getLiquidIsobaricHeatCapacity(T)
            b = [2.49973,0.023454,-3.80136,13.0945,-14.5180];
            Tr = T/SaturatedNitrous.T_CRIT; %Relative temp
            Cp = b(1) * (1 + b(2)*(1-Tr)^-1 + b(3)*(1-Tr) + b(4)*(1-Tr)^2 + b(5)*(1-Tr)^3);
            Cp = Cp*1000; %kJ to J
        end
        
        %Gets the isobaric heat capacity (Cp, in J/KgK) of the saturated vapour at a
        %given temperature (J). Valid up to 30C
        function Cp = getVapourIsobaricHeatCapacity(T)
            b = [132.632,0.052187,-0.364923,-1.20233,0.536141];
            Tr = T/SaturatedNitrous.T_CRIT; %Relative temp
            Cp = b(1) * ( 1+ b(2)*(1-Tr)^(-2/3) + b(3)*(1-Tr)^(-1/3) + b(4)*(1-Tr)^(1/3) + b(5)*(1-Tr)^(2/3) );
            Cp = Cp*1000; %kJ to J
        end
        
        %Gets the isobaric heat capacity (Cp, in J/KgK) of the ideal nitrous oxide gas at a
        %given temperature (J)
        function Cp = getIdealGasIsobaricHeatCapacity(T)
            b = [-0.169903,0.099053,1.20822,-0.248324];
            Tr = T/SaturatedNitrous.T_CRIT; %Relative temp
            Cp = b(1) + b(2)*Tr^(-1/2) + b(3)*Tr^(1/2) + b(4)*Tr;
            Cp = Cp*1000; %kJ to J
        end
        
        %Gets the thermal conductivity (W/mK) of the saturated liquid at a given
        %temperature (K)
        function k = getLiquidThermalConductivity(T)
            b = [72.35,1.5,-3.5,4.5];
            Tr = T/SaturatedNitrous.T_CRIT; %Relative temp
            k = b(1)* (1 + b(2)*(1-Tr)^(1/3) + b(3)*(1-Tr)^(2/3) + b(4)*(1-Tr));
            k = k/1000; %Convert from mW to W
        end
    end
end