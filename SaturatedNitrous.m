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
        T_CRIT = 309.52; %Critical temp of nitrous in kelvin
        P_CRIT = 72.51 * 10^5; %Critical pressure of nitrous, Pa
        DENSITY_CRIT = 452; %Critical density of nitrous, Kg/m^3
    end
    methods (Static) %Static methods: Methods that don't need object instance
        %Get's the specific entropy of saturated nitrous flow mixture of
        %given quality X, temperature T and pressure P
        function s = getSpecificEntropy(X,T,P)
            if(X == 0)
                s = FluidType.NITROUS_LIQUID.getEntropy(T,P);
                return;
            elseif(X == 1)
                s = FluidType.NITROUS_GAS.getEntropy(T,P);
                return;
            end
            sVapour = FluidType.NITROUS_GAS.getEntropy(T,P);
            sLiquid = FluidType.NITROUS_LIQUID.getEntropy(T,P);
            s = sVapour*X + sLiquid*(1-X);
        end
        
        function s = getSaturationSpecificEntropy(X,T)
            P = SaturatedNitrous.getVapourPressure(T);
            s = SaturatedNitrous.getSpecificEntropy(X,T,P);
        end
        
        %Get's the specific enthalpy of saturated nitrous flow mixture of
        %given quality X, temperature T and pressure P
        function h = getSpecificEnthalpy(X,T,P)
            if(X == 0)
               h = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(T,P);
               return;
            elseif(X == 1)
               h = FluidType.NITROUS_GAS.getSpecificEnthalpy(T,P);
               return;
            end
            hVapour = FluidType.NITROUS_GAS.getSpecificEnthalpy(T,P);
            hLiquid = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(T,P);
            h = hVapour*X + hLiquid*(1-X);
        end
        
        %Gets the density of the mixture for a given quality X, temperature
        %T and pressure P
        function rho = getDensity(X,T,P)
              
              rho = NitrousFluidCoolProp.getProperty(FluidProperty.DENSITY,FluidProperty.PRESSURE,P,FluidProperty.VAPOR_QUALITY,X);
%             if(X == 0)
%                 rho = FluidType.NITROUS_LIQUID.getDensity(T,P);
%                 return;
%             elseif(X==1)
%                 rho = FluidType.NITROUS_GAS.getDensity(T,P);
%                 return;
%             end
%             rhoVapour = FluidType.NITROUS_GAS.getDensity(T,P);
%             rhoLiquid = FluidType.NITROUS_LIQUID.getDensity(T,P);
%             rho = 1 / ((1/rhoVapour)*X + (1/rhoLiquid)*(1-X));
        end
        
        %Bubble growth time characteristic as defined in "Modeling feed
        %system flow physics for self pressurising propellants" by Jonny
        %Dyer https://arc.aiaa.org/doi/pdf/10.2514/6.2007-5702. 1 denotes
        %upstream, 2 denotes downstream. Units are SI
        function Tb = getBubbleGrowthTimeCharacteristic(P1,T1,P2)
            vapourPressure = SaturatedNitrous.getVapourPressure(T1);
            rhoLiq = FluidType.NITROUS_LIQUID.getDensity(T1,vapourPressure);
%             disp("T1: "+T1+", rhoLiq: "+rhoLiq+" P2, "+P2+" PVap: "+vapourPressure+" P1: "+P1);
            %vapourPressure = SaturatedNitrous.getVapourPressure(T1);
            Tb = real(sqrt(1.5 * (rhoLiq/(vapourPressure-P2))));
        end
        
        %Residence time (proportional to how long fluid will take to do this) characteristic as defined in "Modeling feed
        %system flow physics for self pressurising propellants" by Jonny
        %Dyer https://arc.aiaa.org/doi/pdf/10.2514/6.2007-5702. 1 denotes
        %upstream, 2 denotes downstream. Units are SI. Length is distance
        %in metres flow is travelling over for this pressure change
        function Tr = getResidenceTime(P1,T1,P2,length)
            vapourPressure = SaturatedNitrous.getVapourPressure(T1);
            rhoLiq = FluidType.NITROUS_LIQUID.getDensity(T1,vapourPressure);
            Tr = real(length*sqrt(rhoLiq/(2*(P1-P2))));
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
            if(isnan(k))
               disp("K IS NAN "+P1+" "+T1+" "+P2); 
               disp(SaturatedNitrous.getBubbleGrowthTimeCharacteristic(P1,T1,P2));
               disp(SaturatedNitrous.getResidenceTime(P1,T1,P2,length));
            end
        end
        
        function [X2,T2,v2,h2,P2] = getDownstreamIncompressiblePressureTemp(X1,P1,T1,G,v1,ignoreTempAndQuality)
            rho = SaturatedNitrous.getDensity(X1,T1,P1); %Upstream density, as incompressible assume is constant
            dP = (0.5*(G.^2 - (rho.^2).*(v1.^2))./rho); %P1 - P2
            P2 = P1 - dP;
            
            v2 = G ./ rho; %G = rho * v
            %Assume total enthalpy constant to figure out T
            h0 = SaturatedNitrous.getSpecificEnthalpy(X1,T1,P1) + (v1.^2)./2;
            %h0 = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(T1,P1) + (v1.^2)./2;
            h2 = h0 - (v2.^2/2); %Downstream enthalpy
            
            if ~exist('ignoreTemp','var') || ~ignoreTempAndQuality
                T2 = SaturatedNitrous.getSaturationTemperature(P2);
                X2 = NitrousFluidCoolProp.getProperty(FluidProperty.VAPOR_QUALITY,FluidProperty.PRESSURE,P2,FluidProperty.SPECIFIC_ENTHALPY,h2);
                if(X2 == -1) %Answer is not on saturation line
                    TSat = T2;
                    T2 = NitrousFluidCoolProp.getProperty(FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P2,FluidProperty.SPECIFIC_ENTHALPY,h2);
                    if(T2 <= TSat)
                        X2 = 0; %Liquid
                    else
                        X2 = 1; %Gas
                    end
                end
            else
                T2 = -1; 
                X2 = 0;
            end
        end
        
        %Get the downstream mass flow rate predicted by an ideal incompressible
        %assumption for a given upstream pressure, temperature and velocity
        %and given downstream pressure. G is mass flow per unit area, all
        %units SI. Temperature given by this also assumes flow is
        %adiabatic and quality given by this assumes isentropic. Set ignoreTempAndQuality to true if you don't want to calculate
        %the temperature and quality downstream (will make this 1000x faster to
        %compute)
        function [X2,T2,v2,h2,G] = getDownstreamIncompressibleMassFlowTemp(X1,P1,T1,P2,v1,ignoreTempAndQuality)
            rho = SaturatedNitrous.getDensity(X1,T1,P1); %Upstream density, as incompressible assume is constant
            G = sqrt(2.*rho.*(P1-P2) + (rho.^2).*v1.^2);
            v2 = G ./ rho; %G = rho * v
            %Assume total enthalpy constant to figure out T
            h0 = SaturatedNitrous.getSpecificEnthalpy(X1,T1,P1) + (v1.^2)./2;
            %h0 = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(T1,P1) + (v1.^2)./2;
            h2 = h0 - (v2.^2/2); %Downstream enthalpy
            
            if ~exist('ignoreTemp','var') || ~ignoreTempAndQuality
                T2 = SaturatedNitrous.getSaturationTemperature(P2);
%                 if(h2 < 399380.7881)
%                     disp("h2 uncorrected was: "+h2);
%                     h2 = 399380.7881+2;
%                 end
                try
                    X2 = NitrousFluidCoolProp.getProperty(FluidProperty.VAPOR_QUALITY,FluidProperty.PRESSURE,P2,FluidProperty.SPECIFIC_ENTHALPY,h2);
                catch
                    X2 = NitrousFluidCoolProp.getProperty(FluidProperty.VAPOR_QUALITY,FluidProperty.PRESSURE,P2,FluidProperty.SPECIFIC_ENTHALPY,h2+2);
                end
                
                if(X2 == -1) %Answer is not on saturation line
                    TSat = T2;
                    try
                        T2 = NitrousFluidCoolProp.getProperty(FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P2,FluidProperty.SPECIFIC_ENTHALPY,h2);
                    catch
                        T2 = NitrousFluidCoolProp.getProperty(FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P2,FluidProperty.SPECIFIC_ENTHALPY,h2+2);
                    end
                    if(T2 <= TSat)
                        X2 = 0; %Liquid
                    else
                        X2 = 1; %Gas
                    end
                end
            else
                T2 = -1; 
                X2 = 0;
            end
%             function err = T2Err(TGuess) %Error in guess for downstream temperature (to match the enthalpy)
%                 h = SaturatedNitrous.getSpecificEnthalpy(X1,real(TGuess),P2);
%                 %h =  FluidType.NITROUS_LIQUID.getSpecificEnthalpy(real(TGuess),P2);
%                 err = h-h2; %0 if downstream enthalpy matches
%             end
        end
        
        function [TSat,PToDropToForSaturationLine] = getWhereIsentropicIntersectWithSaturationCurve(T,P)
            TSatUpstream = SaturatedNitrous.getSaturationTemperature(P);
            liquidUpstream = T<TSatUpstream;
            if(abs(1-T/TSatUpstream) < 0.0001)
               TSat = T;
               PToDropToForSaturationLine = P;
               return;
            end
            XUpstream = 1;
            if(liquidUpstream)
               XUpstream = 0;
               s1 = FluidType.NITROUS_LIQUID.getEntropy(T,P);
            else
               s1 = FluidType.NITROUS_GAS.getEntropy(T,P);
            end
            
            errSatT = @(T) s1 - SaturatedNitrous.getSaturationSpecificEntropy(XUpstream,real(T));
            TSat = real(betterfzero(errSatT,T,183,SaturatedNitrous.T_CRIT,1e-3));
            PToDropToForSaturationLine = SaturatedNitrous.getVapourPressure(TSat);
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
                %h1 = SaturatedNitrous.getSpecificEnthalpy(X1,T1,P1); %Enthalpy upstream
                [T1,vRealFlow,~] = RealFlow.getIsentropicTempVelocity(P1,T1,v1,vapourPressure,FluidType.NITROUS_LIQUID,1);
                %disp("vRealFlow: "+vRealFlow);
                %[T1,~,~] = RealFlow.getIsentropicTempVelocity(P1,T1,T1,FluidType.NITROUS_LIQUID.getCp(T1,P1),vapourPressure,FluidType.NITROUS_LIQUID);
%                 h2 = SaturatedNitrous.getSpecificEnthalpy(X1,T1,vapourPressure); %Enthalpy after pressure drop
%                 KEGain = h1-h2; %specific KE gain is difference in fluid enthalpies upstream and downstream
%                 v1 = sqrt(2*(0.5*v1.^2 + KEGain)); %Energy balance, assuming no other losses
                v1 = real(vRealFlow); 
                P1 = vapourPressure;
                if P1 == P2 %Flow never reached saturation line
                    h2 = SaturatedNitrous.getSpecificEnthalpy(X1,T1,vapourPressure); %Enthalpy after pressure drop
                    X2 = X1;
                    T2 = T1;
                    if(imag(T2)~=0)
                       error('T2 became complex (1)'); 
                    end
                    v2 = v1;
                    %h2 = h2;
                    rhoLiq = FluidType.NITROUS_LIQUID.getDensity(T2,P2);
                    G = rhoLiq * v2;
                    return;
                end
            end      
            
            %While a vapour not yet at saturation line (P < Vapour
            %pressure within 4%) use RealFlow to handle transport
            if(X1>0.99 && ((P1-SaturatedNitrous.getVapourPressure(T1))/SaturatedNitrous.getVapourPressure(T1)) < -0.04)
                %disp("VAPOUR ONLY: "+((P1-SaturatedNitrous.getVapourPressure(T1))/SaturatedNitrous.getVapourPressure(T1)));
                s1 = FluidType.NITROUS_GAS.getEntropy(T1,P1);
                %Find what pressure need to drop to so that intersect with
                %saturation line
                %distFromSaturationLine = @(P) FluidType.NITROUS_GENERAL.getTemperatureFromPressureEntropy(real(P),s1) - SaturatedNitrous.getSaturationTemperature(real(P));
                %dDdP = (distFromSaturationLine(P1)-distFromSaturationLine(P1-100)) / (-100);
                %approxdP = (1/dDdP) * distFromSaturationLine(P1);
                %PToDropToForSaturationLine = P1+approxdP;
                errSatT = @(T) s1 - SaturatedNitrous.getSaturationSpecificEntropy(1,real(T));
                try
                    TSat = real(betterfzero(errSatT,T1,183,SaturatedNitrous.T_CRIT,1e-3));
                catch err
                    %Search failed, assume TSat is 183K. Eg. probably never
                    %going to reach it
                    TSat = 183;
                end
                PToDropToForSaturationLine = SaturatedNitrous.getVapourPressure(TSat);
%                 PToDropToForSaturationLine = NitrousFluidCoolProp.getProperty(FluidProperty.TEMPERATURE,'Smass|twophase',s1,FluidProperty.VAPOR_QUALITY,1);
              %   disp("Saturation P: "+PToDropToForSaturationLine);
                if PToDropToForSaturationLine < P2 %Fluid remains as a gas, no need to use NHNE model
                    [T2,v2,~] = RealFlow.getIsentropicTempVelocity(P1,T1,v1,P2,FluidType.NITROUS_GAS,1);
                    X2 = X1;
                    h2 = FluidType.NITROUS_GAS.getEnthalpy(T2,P2);
                    G = v2 * FluidType.NITROUS_GAS.getDensity(T2,P2); %G = rho * v
                    if(imag(T2)~=0)
                       error('T2 became complex (2)'); 
                    end
                    return;
                end
                %Drop fluid to saturation line then continue with NHNE
                %model
                [T1,v1,~] = RealFlow.getIsentropicTempVelocity(P1,T1,v1,PToDropToForSaturationLine,FluidType.NITROUS_GAS,1);
                P1 = PToDropToForSaturationLine;
            end
            
            k = SaturatedNitrous.getNonEquilibriumFlowParameter(P1,T1,P2,length);
            
            incompressibleCoeff = (1-1/(1+k)); %From correction to dyer by solomon
            hemCoeff = (1/(1+k)); %From correction to dyer by solomon
            
            if(imag(T1)~=0)
                error('T1 became complex (1)');
            end
            
            %Calculate flow both incompressibly and with HEM
            [X2Inc,T2Inc,v2Inc,h2Inc,GInc] = SaturatedNitrous.getDownstreamIncompressibleMassFlowTemp(X1,P1,T1,P2,v1);
            [X2Hem,T2Hem,v2Hem,h2Hem,GHem] = SaturatedNitrous.getDownstreamIsentropicSaturatedHEMFlowCond(X1,T1,P1,P2,v1);
            %Weighting of flow properties as in Dyer
            X2 = incompressibleCoeff * X2Inc + hemCoeff * X2Hem;
            T2 = incompressibleCoeff * T2Inc + hemCoeff * T2Hem;
            %disp("x2: "+X2+" T2: "+T2+" hem X2: "+X2Hem+" inc X2: "+X2Inc+" T2Hem: "+T2Hem+" T2Inc: "+T2Inc);
            %v2 = incompressibleCoeff * v2Inc + hemCoeff * v2Hem;
            h2 = incompressibleCoeff * h2Inc + hemCoeff * h2Hem;
            G = incompressibleCoeff * GInc + hemCoeff * GHem;
%             disp("P2: "+P2+", GInc: "+GInc+" GHem: "+GHem);
            %G = rho*v

            if(X2 < 0 && X2 > -0.000001)
               X2 = 0; 
            elseif(X2 > 1 && X2 < 1.000001)
                   X2 = 1; 
            end
            
            rhoDownstream = SaturatedNitrous.getDensity(X2,T2,P2);
            v2 = G ./ rhoDownstream;
            if(imag(T2)~=0)
                disp(k);
                disp(P1);
                disp(T1);
                disp(P2);
                disp(length);
                error('T2 became complex (3)');
            end
        end
        
        function [X2,T2,v2,h2,P2] = getDownstreamIsentropicSaturationHEMFlowPressure(X1,T1,P1,G,v1)
            rho1 = SaturatedNitrous.getDensity(X1,T1,P1);
            P2Old = -1;
            rho2 = rho1; %Initially assume downstream density is the same, refine this through iteration
            X2 = X1; %Initial value, will be changed
            T2 = T1; %Initial value will be changed
            s1 = SaturatedNitrous.getSpecificEntropy(X1,T1,P1); %Entropy upstream
            h1 = SaturatedNitrous.getSpecificEnthalpy(X1,T1,P1); %Enthalpy upstream
            iter = 0;
            P2 = P1; %Initial value, will change
            while(abs(P2-P2Old) > 0.5 && iter < 100) %Until converge or 100 iterations
                P2Old = P2;
                iter = iter+1;
                v2 = G ./ rho2;
                KEGain = 0.5.*v2.^2 - 0.5.*v1.^2;
                h2 = h1-KEGain;
                
                hErr = @(P) h2 - NitrousFluidCoolProp.getProperty(FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.PRESSURE,P,FluidProperty.DENSITY,rho2);
%                 tic;
                %P2 = NitrousFluid.getPressureFromDensityEnthalpy(rho2,h2);
                P2 = fzero(hErr,P1);
%                 toc;
%                 tic;
%                 P22 = NitrousFluid.getPressureFromDensityEnthalpy(rho2,h2);
%                 
%                 P2 = fzero(hErr,P22);
%                 toc;
%                 disp(P2+" "+P22+" "+abs(P22-P2));
                 
                %P2 = NitrousFluidCoolProp.getProperty(FluidProperty.PRESSURE,FluidProperty.SPECIFIC_ENTHALPY,h2,FluidProperty.TEMPERATURE,T2);
                T2 = SaturatedNitrous.getSaturationTemperature(P2); %Drops to saturation temp of new pressure
                X2 = NitrousFluidCoolProp.getProperty(FluidProperty.VAPOR_QUALITY,FluidProperty.PRESSURE,P2,FluidProperty.SPECIFIC_ENTROPY,s1);
                if(X2 == -1) %NOT Saturated
                    TSat = T2;
                    T2 = NitrousFluidCoolProp.getProperty(FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P2,FluidProperty.SPECIFIC_ENTROPY,s1);
                    if(T2 <= TSat)
                        X2 = 0; %Liquid
                    else
                        X2 = 1; %Gas
                    end
                    rho2 = NitrousFluidCoolProp.getProperty(FluidProperty.DENSITY,FluidProperty.TEMPERATURE,T2,FluidProperty.PRESSURE,P2);
                else
                    rho2 = NitrousFluidCoolProp.getProperty(FluidProperty.DENSITY,FluidProperty.PRESSURE,P2,FluidProperty.VAPOR_QUALITY,X2);
                end
            end
%             disp("Iter num: "+iter);
%             [X2Test,T2Test,v2Test,h2Test,GTest] = SaturatedNitrous.getDownstreamIsentropicSaturatedHEMFlowCond(X1,T1,P1,P2,v1);
%             disp(GTest);
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
                [T1,~,~] = RealFlow.getIsentropicTempVelocity(P1,T1,v1,P2,FluidType.NITROUS_GENERAL);
                %[T1,~,~] = RealFlow.getIsentropicTempVelocity(P1,T1,T1,FluidType.NITROUS_LIQUID.getCp(T1,P1),vapourPressure,FluidType.NITROUS_LIQUID);
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
            X2 = NitrousFluidCoolProp.getProperty(FluidProperty.VAPOR_QUALITY,FluidProperty.PRESSURE,P2,FluidProperty.SPECIFIC_ENTROPY,s1);
            %X2 = real(fzero(@(X) real(SaturatedNitrous.getSpecificEntropy(real(X),T2,P2)-s1),X1)); %Find quality downstream that satisfies isentropic condition
            if X2 < 0 || X2 > 1
                if(X2 == -1)
                   %Not saturated
                   TSat = T2;
                   T2 = NitrousFluidCoolProp.getProperty(FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P2,FluidProperty.SPECIFIC_ENTROPY,s1);
                   if(T2 <= TSat)
                       X2 = 0; %Liquid
                   else
                       X2 = 1; %Gas
                   end
                elseif(X2 >1 && X2<1.000001)
                    X2 = 1; %Gas 
                elseif(X2 <0 && X2>-0.000001)
                    X2 = 0; %Liquid
                else
                    disp("INVALID X2 of "+X2);
                    error('INVALID x2 encountered'); 
                end
%                 disp("X2: "+X2);
%                 disp("P1: "+P1);
%                 disp("P2: "+P2);
%                 disp("X1: "+X1);
            end
            h2 = SaturatedNitrous.getSpecificEnthalpy(X2,T2,P2); %Gets the specific enthalpy downstream
            KEGain = h1-h2; %specific KE gain is difference in fluid enthalpies upstream and downstream
            v2 = sqrt(2*(0.5*v1.^2 + KEGain)); %Energy balance, assuming no other losses
            rhoMixture = SaturatedNitrous.getDensity(X2,T2,P2);
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
            
            T = NitrousFluidCoolProp.getProperty(FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P,FluidProperty.VAPOR_QUALITY,0);
%             T = fzero(@(T) (real(SaturatedNitrous.getVapourPressure(real(T)) - P)),290);
%             disp(T0+" "+T);
%             if abs(imag(SaturatedNitrous.getVapourPressure(real(T)))) > 0
%                 warning('No saturation temperature for given pressure, fluid cannot be liquid'); 
%                T = -1; %Does not exist 
%             end
        end
        
        function P = getVapourPressureFromFormula(T)
            TRel = T / SaturatedNitrous.T_CRIT; %Relative temp
            b = [-6.71893,1.35966,-1.3779,-4.051]; %Constants in equation
            %Equation from paper
            logPressureRatio = (1/TRel) * ( b(1)*(1-TRel) + b(2)*(1-TRel)^1.5 + b(3)*(1-TRel)^2.5 + b(4)*(1-TRel)^5 );
            P = exp(logPressureRatio) * SaturatedNitrous.P_CRIT; %Calculate vapour pressure
        end
        
        %Gets the vapour pressure (Pa) of saturated nitrous at a given
        %temperature (K)
        function P = getVapourPressure(T)
            %Polynomial fit to coolprop's output for speed enhancement
            polyCoeffs = [6.30913881904569e-14,-1.51190519529204e-10,1.62556872417980e-07,-0.000103259603725827,0.0429133163140751,-12.1912294858819,2397.61266586979,-322316.427676396,28344892.9083561,-1472430566.16391,34309300122.4424];
            P = polyval(polyCoeffs,T);
%             P = NitrousFluidCoolProp.getProperty(FluidProperty.PRESSURE,FluidProperty.TEMPERATURE,T,FluidProperty.VAPOR_QUALITY,0);
%             TRel = T / SaturatedNitrous.T_CRIT; %Relative temp
%             b = [-6.71893,1.35966,-1.3779,-4.051]; %Constants in equation
%             %Equation from paper
%             logPressureRatio = (1/TRel) * ( b(1)*(1-TRel) + b(2)*(1-TRel)^1.5 + b(3)*(1-TRel)^2.5 + b(4)*(1-TRel)^5 );
%             P = exp(logPressureRatio) * SaturatedNitrous.P_CRIT; %Calculate vapour pressure
        end
        
        %Gets the liquid density (Kg/m^3) of saturated nitrous at a given
        %temperature (K)
        function rho = getLiquidDensity(T)
%             rho = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.DENSITY, FluidProperty.TEMPERATURE,T,FluidProperty.VAPOR_QUALITY,0);
            TRel = T / SaturatedNitrous.T_CRIT; %Relative temp
            b = [1.72328,-0.83950,0.51060,-0.10412]; %Constants in equation
            %Equation from paper
            logDensityRatio =  b(1)*(1-TRel)^(1/3) + b(2)*(1-TRel)^(2/3) + b(3)*(1-TRel)^(1) + b(4)*(1-TRel)^(4/3);
            rho = exp(logDensityRatio) * SaturatedNitrous.DENSITY_CRIT; %Calculate vapour pressure
        end
        
        %Gets the vapour density (Kg/m^3) of saturated nitrous at a given
        %temperature (K)
        function rho = getVapourDensity(T)
%               rho = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.DENSITY, FluidProperty.TEMPERATURE,T,FluidProperty.VAPOR_QUALITY,1);
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