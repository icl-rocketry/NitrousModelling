%Prefer to use NitrousFluidCoolProp.x directly
%Internally delegates to coolprop in some places, and in other places
%attempts to interpolate between NIST data. Be carefuld some NIST and
%coolprop values are not directly compatible (Eg. energy values)
%
%Class to get properties of nitrous gas and liquid (N2O).
%Depends on folder 'nitrousRawData' and files contained within
%Data from https://wtt-pro.nist.gov to determine values
%By Eddie Brown
classdef NitrousFluid
    properties(Constant)
        SUPPRESS_OUT_OF_RANGE_WARNINGS = true;
    end
    methods (Static,Access=private) %Static private methods: Methods that don't need object instance that can't be called remotely
        
        %Function to linearly interpolate fastly with one column acting as
        %independent variable (Eg. input, eg. temp) and one column acting
        %as dependent variable (Eg. output, eg. Enthalpy)
        function val = oneColInterp(data,indepCol1,depCol,indepVal1)
            if (length(indepVal1) > 1) %Handle vector input
                val = zeros(1,length(indepVal1));
                for j=1:length(indepVal1)
                   val(j) = NitrousFluid.oneColInterp(data,indepCol1,depCol,indepVal1(j)); 
                end
                return;
            end
            %Find indexes of positions in list either side of correct value
            %for indepVal1
            [~,indexLowerBound,indexUpperBound,~] = NitrousFluid.binarySearch(data,indepCol1,indepVal1);
            lowerBoundIndepVal = data(indexLowerBound,indepCol1); %Lower bound of indep val
            upperBoundIndepVal = data(indexUpperBound,indepCol1);
            lowerBoundDepVal = data(indexLowerBound,depCol); %Lower bound of dependent val
            upperBoundDepVal = data(indexUpperBound,depCol);
            
            %Come up with a third index to do a cubic interpolation
            thirdIndex = indexLowerBound -1; %Data point before bounded region
            if(thirdIndex < 1) %If doesn't exist use data point after bounded region
               thirdIndex = indexUpperBound+1; 
            end
            dataSize = size(data);
            if(thirdIndex > dataSize(1)) %No possible third data point
                thirdIndex = indexLowerBound; %Make third data point the same as the first
            end
            thirdPtIndepVal = data(thirdIndex,indepCol1);
            thirdPtDepVal = data(thirdIndex, depCol);
            
            %Cubic interpolate
%             if indepVal1 < lowerBoundIndepVal || indepVal1 > upperBoundIndepVal
%                disp("Interpolating outside of range"); 
%                disp("1: "+lowerBoundIndepVal+" > "+lowerBoundDepVal);
%                disp("2: "+upperBoundIndepVal+" > "+upperBoundDepVal);
%                disp("3: "+thirdPtIndepVal+" > "+thirdPtDepVal);
%                disp("Indep val: "+indepVal1);
% %                val = NitrousFluid.linearInterp(lowerBoundIndepVal,upperBoundIndepVal,lowerBoundDepVal,upperBoundDepVal,indepVal1);
% %                 return;
%             end
            %val = NitrousFluid.linearInterp(lowerBoundIndepVal,upperBoundIndepVal,lowerBoundDepVal,upperBoundDepVal,indepVal1);
            val = NitrousFluid.cubicInterp(lowerBoundIndepVal,upperBoundIndepVal,thirdPtIndepVal,lowerBoundDepVal,upperBoundDepVal,thirdPtDepVal, indepVal1);
        end
        
        %Function to linearly interpolate fastly with two columns acting as
        %independent variables (eg. inputs, eg. x) and one as the dependent
        %variable (eg. output, eg. y).
        %Eg. val = twoColLinearInterp(data,indexOfColOfTemp,indexOfColOfPressure,indexOfColOfValueYouWant,temp,pressure);
        %ASSUMES that the data is ordered by indepCol1 and that rows with
        %the same value for indepCol1 are ordered by indepCol2
        function val = twoColInterp(data,indepCol1,indepCol2,depCol,indepVal1,indepVal2)
            if (length(indepVal1) > 1) %Handle vector input
                if(length(indepVal1) ~= length(indepVal2))
                    error('Two col linear interp is not capable of handling input given as differently sized vectors');
                end
                val = zeros(1,length(indepVal1));
                for j=1:length(indepVal1)
                   val(j) = NitrousFluid.twoColInterp(data,indepCol1,indepCol2,depCol,indepVal1(j),indepVal2(j)); 
                end
                return;
            end
            %Find indexes of positions in list either side of correct value
            %for indepVal1
            [indep1StartIndexLower,indep1EndIndexLower,indep1StartIndexUpper,indep1EndIndexUpper] = NitrousFluid.binarySearch(data,indepCol1,indepVal1);
            indep1Lower = data(indep1EndIndexLower,indepCol1); %The value of the lower bound for the first independent variable. Eg. if var=T and T we want is 23.5, this might be 23
            indep1Upper = data(indep1StartIndexUpper,indepCol1); %The value of the upper bound for the first independent variable. Eg. if var=T and T we want is 23.5, this might be 24
            
            dataAtIndep1Lower = data(indep1StartIndexLower:indep1EndIndexLower,:); %All the data for when indepdent column 1 is the lower bound value. Eg. this is all values where T=T_lower_bound
            dataAtIndep1Upper = data(indep1StartIndexUpper:indep1EndIndexUpper,:); %All the data for when indepdent column 1 is the upper bound value. Eg. this is all values where T=T_upper_bound
            
            %Interpolated second column values at each of the first columns
            %required for the second interpolation (Use these values for
            %2nd interpolation)
            valLower = NitrousFluid.oneColInterp(dataAtIndep1Lower,indepCol2,depCol,indepVal2);
            valUpper = NitrousFluid.oneColInterp(dataAtIndep1Upper,indepCol2,depCol,indepVal2);
%             disp(data(:,2:3));
            
            %Need to find a third indep1 for cubic interpolation
            dataSize = size(data);
            if(indep1StartIndexLower - 1 < 1) %If no data points for earlier values of indep1
                if(indep1EndIndexUpper+1>dataSize(1)) %No data points for later values of indep1
                    %No other data points to grab, so just linearly
                    %interpolate
                    %Linearly interpolate on col1 to get final interpolated result
                    val = NitrousFluid.linearInterp(indep1Lower,indep1Upper,valLower,valUpper,indepVal1);
                    return;
                end
                
                %Use data points for later values of indep1
                thirdPtStartIndex = indep1EndIndexUpper+1;
                thirdPtEndIndex = thirdPtStartIndex; %Init end index of range for data which has same value of indep1 at start index
                while (thirdPtEndIndex <= dataSize(1) && data(thirdPtEndIndex,indepCol1) == data(thirdPtStartIndex,indepCol1)) %Go through data until find different value of indep1
                    thirdPtEndIndex = thirdPtEndIndex+1;
                end
                thirdPtEndIndex = thirdPtEndIndex-1; %Make the end index the final index that had the same value of indep1
            else %Data points for earlier values of indep1, use these for our third interpolation point
                thirdPtEndIndex = indep1StartIndexLower - 1;
                thirdPtStartIndex = thirdPtEndIndex; %Init start index of range for data which same value of indep1 at the end index
                while (thirdPtStartIndex >= 1 && data(thirdPtStartIndex,indepCol1) == data(thirdPtEndIndex,indepCol1)) %Go through data until find different value of indep1
                    thirdPtStartIndex = thirdPtStartIndex - 1;
                end
                thirdPtStartIndex = thirdPtStartIndex+1;
            end
            dataAtThirdIndep1 = data(thirdPtStartIndex:thirdPtEndIndex,:);
            thirdIndep1 = data(thirdPtStartIndex,indepCol1);
            valThird = NitrousFluid.oneColInterp(dataAtThirdIndep1,indepCol2,depCol,indepVal2);
            
            %Interpolate on col1 to get final interpolated result
            val = NitrousFluid.cubicInterp(indep1Lower,indep1Upper,thirdIndep1,valLower,valUpper,valThird,indepVal1);
            %val = NitrousFluid.linearInterp(indep1Lower,indep1Upper,valLower,valUpper,indepVal1);
        end
        
        %Simple linear interpolation function
        function val = linearInterp(x1,x2,y1,y2,x)
            if x2-x1 == 0 %If only one data point
                val = y1; %Return value of single data point
                return;
            end
            val = ((x-x1)./(x2-x1)).*(y2-y1) + y1;
        end
        
        %Simple cubic interpolation function
        function val = cubicInterp(x1,x2,x3,y1,y2,y3,x)
            if (x3 == x1 || x3 == x2) %No third data point, so are re-using a 2nd data point
                val = NitrousFluid.linearInterp(x1,x2,y1,y2,x); %Linearly interpolate instead
                return;
            end
            xVals = [x1,x2,x3];
            yVals = [y1,y2,y3];
            val = spline(xVals,yVals,x); %Piecewise Cubic Interpolating Polynomial
        end
        
        %Function to use modified binary search on column with index 'colIndex' within ordered dataset 'data'
        %for the indexes of the values closest to 'val' (either side of it). Returns index where first
        %appears and index where last appears for each of these values.
        %Purpose is for then linearly interpolating
        function [startIndexLower,endIndexLower,startIndexUpper,endIndexUpper] = binarySearch(data,colIndex,val)
            [numRows,~] = size(data);
            domainStart = 1; %Index of the start of the list used for binary search
            domainEnd = numRows; %Index of the end of the list used for binary search
            increasing = true;
            if (length(data) > 1 && data(2,colIndex) < data(1,colIndex)) %If data in descending order
               increasing = false; 
            end
            found = false;
            while ~found
                midIndex = ceil((domainStart + domainEnd) / 2);
                
                %Check if 'found' at midIndex. In this weird case 'found'
                %means value found <= 'val', but next largest element is >
                %val - so that the value being found by the search is the lower of the
                %two required for linearly interpolating
                valueFound = data(midIndex,colIndex);
                if (increasing && valueFound <= val) || (~increasing && valueFound >= val)
                    %check if found
                    for i=midIndex+1:numRows %Iterate through list starting at midIndex
                        if data(i,colIndex) == valueFound && i~=numRows %If value is the same as at midIndex and we are not yet at the end of the list
                            continue; %Continue iterating over the list
                        end
                        %i is now the index of the first element in the
                        %list after midIndex that is larger, or we have
                        %reached the end of the list
                        if (data(i,colIndex) >= val && increasing) || (~increasing && data(i,colIndex) <= val) %val is between value at position i and value at position (i-1), yay!
                            found = true;
                            startIndexUpper = i; %We know that i is the index of the left-most smallest value larger than val in the list
                            endIndexLower = i-1; %We know that i-1 is the index of the right-most largest value smaller than val in the list
                            break; %Break the for loop
                        end
                    end
                    if ~found %If was not found then valueFound is too small
                        domainStart = midIndex + 1;
                        if domainStart > domainEnd %Search failed to find value that matched criteria
                            break; %End the search
                        end
                    end
                else %valueFound is > val
                    domainEnd = midIndex - 1; %Move domain to exclude segment of list too large
                    if domainEnd < domainStart %Search failed to find value that matched criteria
                        break; %End the search
                    end
                end
            end
            
            if ~found %If were unable to find endIndexLower,startIndexUpper that are satisfactory
                MAX_INTERP_OUTSIDEOFRANGE_NO_WARN = 1000; %Maximum amount to interpret by outside of the data range before warning
                if numRows >= 4 %If  rows or more
                    %Maximum interp without warn is 5*spacing between 2nd and
                    %3rd values in list
                    MAX_INTERP_OUTSIDEOFRANGE_NO_WARN  = 5*abs(data(3,colIndex) - data(2,colIndex));
                end
                if (increasing && val >= data(numRows,colIndex)) || (~increasing && val <= data(numRows,colIndex)) %Val is bigger than largest in dataset
                    amountOutsideRange = abs(val - data(numRows,colIndex));
                    if ~NitrousFluid.SUPPRESS_OUT_OF_RANGE_WARNINGS && amountOutsideRange > MAX_INTERP_OUTSIDEOFRANGE_NO_WARN
                        warning(['Value ',num2str(val),' is outside of dataset, closest value is ',num2str(data(numRows,colIndex)),'. Interpolation may be less accurate']);
                    end
                    endIndexUpper = numRows;
                    %Find startIndexUpper
                    for startIndexUpper=endIndexUpper:-1:1 %Iterate from endIndexUpper to start of list
                        %Once reach a point where the value of the position
                        %before
                        %startIndexUpper is a different value, or end of list
                        if startIndexUpper-1 < 1
                            endIndexLower = startIndexUpper;
                            break; %Stop iterating, now startIndexUpper is at correct value
                        end
                        if startIndexUpper-2 < 1 || data(startIndexUpper-1,colIndex) ~= data(startIndexUpper,colIndex)
                            endIndexLower = startIndexUpper-1;
                            break; %Stop iterating, now startIndexUpper is at correct value
                        end
                    end
                elseif (increasing && val <= data(1,colIndex)) || (~increasing && val >= data(1,colIndex)) %Val is smaller than smallest in dataset
                    amountOutsideRange = abs(data(1,colIndex)-val);
                    if ~NitrousFluid.SUPPRESS_OUT_OF_RANGE_WARNINGS && amountOutsideRange > MAX_INTERP_OUTSIDEOFRANGE_NO_WARN
                        warning(['Value ',num2str(val),' is outside of dataset, closest value is ',num2str(data(1,colIndex)),'. Interpolation may be less accurate']);
                    end
                    startIndexLower = 1;
                    %Find endIndexLower
                    for endIndexLower=startIndexLower:1:numRows %Iterate from startIndexLower to endof list
                        %Once reach a point where the value at the index of the position after
                        %endIndexUpper is a different value, or end of list
                        if endIndexLower+1 > numRows
                            startIndexUpper = endIndexLower;
                            break; %Stop iterating, now endIndexLower is at correct value
                        end
                        if endIndexLower+2 > numRows || data(endIndexLower+1,colIndex) ~= data(endIndexLower,colIndex)
                            startIndexUpper = endIndexLower+1;
                            break; %Stop iterating, now endIndexLower is at correct value
                        end
                    end
                else
                    disp(data);
                   error(['For some reason unable to find value ',num2str(val),' within dataset where min is ',num2str(data(1,colIndex)),' and max is ',num2str(data(numRows,colIndex))]); 
                end
               %error(['Unable to find values for interpolating between for ',num2str(val),' within column indexed ',num2str(colIndex),'! Dataset probably does not contain this value!']);
            end
            
            %As execution reached here, found=true
            %Find startIndexLower
            for startIndexLower=endIndexLower:-1:1 %Iterate from endIndexLower to 1 backwards
                %Once reach a point where the value at the index of the position before
                %startIndexLower is a different value, or start of list
                if startIndexLower-1 == 0 || data(startIndexLower-1,colIndex) ~= data(startIndexLower,colIndex)
                   break; %Stop iterating, now startIndexLower is at correct value
                end
            end
            
            %Find endIndexUpper
            for endIndexUpper=startIndexUpper:1:numRows %Iterate from startIndexUpper to endof list
                %Once reach a point where the value at the index of the position after
                %endIndexUpper is a different value, or end of list
                if endIndexUpper+1 > numRows || data(endIndexUpper+1,colIndex) ~= data(endIndexUpper,colIndex)
                   break; %Stop iterating, now endIndexUpper is at correct value
                end
            end
        end
        
        %Function to read contents of file, ignoring the first 2 lines, that contains 4 columns and then put into matrix and then
        %return. Implements a cache internally for much faster execution
        %speeds when repeated calls are necessary
        function data = getDataFromFile(fName,numCols)
           persistent cachedData; %Map that caches contents read from files, to speed up execution massively
           if isempty(cachedData) %If variable not initialized
               cachedData = containers.Map('KeyType','char','ValueType','any'); %Initialize it to our map
           end
           if ~exist('numCols','var')
               numCols = 4; 
           end
           if ~isKey(cachedData,fName) %If map does not contain the contents of the file
               if(endsWith(fName,'.mat'))
                   fileData = load(fName);
                   cachedData(fName) = fileData.data;
                   data = fileData.data;
                   return;
               end
               
               %Load data from the file and put into map
               fileHandle = fopen(fName,'r'); %Open file with read perms
               
               fileData = zeros(1,numCols);
               fgetl(fileHandle); %Ignore first line
               fgetl(fileHandle); %Ignore second line
               lineNum = 1;
               %Regex (.+?)\s(.+?)\s(.+?)\s(.+?) for parsing columns,
               %lookup regular expressions in programming if not sure
               %what this is
               regex = '(.+?)';
               for j=1:numCols-1
                   regex = [regex,'\s+(.+'];
                   if(j ~= numCols-1)
                       regex = [regex,'?)'];
                   else
                       regex = [regex,')'];
                   end
               end
               while true %Until loop breaks
                   lineRead = fgetl(fileHandle); %Read next line from file
                   if ~ischar(lineRead) %If this line does not exist (reached end of file)
                       break; %Exit loop
                   end
                   
                   tokens = regexp(lineRead,regex,'tokens'); %Capture groups from the regex as an array
                   for i=1:numCols %For each col, 1 to 4
                        %Put into matrix the numeric value captured by this
                        %token
                        fileData(lineNum,i) = str2double(tokens{1}{i});
%                         if isnan(fileData(lineNum,i))
%                            disp(tokens{1}{i});
%                         end
                   end
                   lineNum = lineNum+1;
               end
               fclose(fileHandle);
               cachedData(fName) = fileData;
           end
           data = cachedData(fName); %Data is what we loaded from this file
        end
    end
    methods (Static) %Static methods: Methods that don't need object instance
        
        %Function to get the adiabatic compressibility (1/Pa) for the gas
        %at a given Temp (K) and Pressure (Pa)
        function val = getGasAdiabaticCompressibility(T,P)
            P1 = P/1000; %Need gas in kPa using tabulated data
            data = NitrousFluid.getDataFromFile(['nitrousRawData',filesep,'gasAdiabaticCompressibility.txt']); 
            val = NitrousFluid.twoColInterp(data,1,2,3,T,P1); %Returns isothermal compressibility in 1/kPa
            val = val/1000; %Convert to 1/Pa
        end
        
        %Function to get the adiabatic compressibility (1/Pa) for the
        %liquid
        %at a given Temp (K) and Pressure (Pa)
        function val = getLiquidAdiabaticCompressibility(T,P)
            P1 = P/1000; %Need gas in kPa using tabulated data
            data = NitrousFluid.getDataFromFile(['nitrousRawData',filesep,'liquidAdiabaticCompressibility.txt']); 
            val = NitrousFluid.twoColInterp(data,1,2,3,T,P1); %Returns isothermal compressibility in 1/kPa
            val = val/1000; %Convert to 1/Pa
        end
        
        %Function to get the isothermal compressibility (1/Pa) for the
        %gas
        %at a given Temp (K) and Pressure (Pa)
        function val = getGasIsothermalCompressibility(T,P)
            P1 = P/1000; %Need gas in kPa using tabulated data
            data = NitrousFluid.getDataFromFile(['nitrousRawData',filesep,'gasIsothermalCompressibility.txt']); 
            val = NitrousFluid.twoColInterp(data,1,2,3,T,P1); %Returns isothermal compressibility in 1/kPa
            val = val/1000; %Convert to 1/Pa
        end
        
        %Function to get the isothermal compressibility (1/Pa) for the
        %liquid
        %at a given Temp (K) and Pressure (Pa)
        function val = getLiquidIsothermalCompressibility(T,P)
            P1 = P/1000; %Need gas in kPa using tabulated data
            data = NitrousFluid.getDataFromFile(['nitrousRawData',filesep,'liquidIsothermalCompressibility.txt']); 
            val = NitrousFluid.twoColInterp(data,1,2,3,T,P1); %Returns isothermal compressibility in 1/kPa
            val = val/1000; %Convert to 1/Pa
        end
        
        %Function to get the molar mass of nitrous oxide in Kg/mol
        function val = getMolarMass()
           val = 44.013 / 1000.0; %Molar mass of 44.013 g/mol
        end
        
        %Function to get the gas constant of nitrous oxide in SI units (J/KgK)
        function val = getGasConstant()
           val = 8.314472 / NitrousFluid.getMolarMass(); %Universal gas constant divide by molar mass
        end
        
        %Function to get the speed of sound of the nitrous oxide in m/s for
        %a given temperature (K) and pressure (P), automatically assuming
        %the correct state
        function c = getSpeedOfSound(T,P)
           c = NitrousFluidCoolProp.getProperty(FluidProperty.SPEED_OF_SOUND,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        %Function to get the speed of sound of the nitrous oxide gas in m/s
        %for a given temperature (K) and pressure (P)
        function c = getGasSpeedOfSound(T,P)
           %Use c^2 = 1 / (rho * isentropic compressibility);
           c = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.SPEED_OF_SOUND,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
%            betaS = NitrousFluid.getGasAdiabaticCompressibility(T,P);
%            rho = NitrousFluid.getGasDensity(T,P);
%            c = sqrt(1 / (rho * betaS));
        end
        
        %Function to get the speed of sound of the nitrous oxide liquid in m/s
        %for a given temperature (K) and pressure (P)
        function c = getLiquidSpeedOfSound(T,P)
            c = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPEED_OF_SOUND,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
           %Use c^2 = 1 / (rho * isentropic compressibility);
%            betaS = NitrousFluid.getLiquidAdiabaticCompressibility(T,P);
%            rho = NitrousFluid.getLiquidDensity(T,P);
%            c = sqrt(1 / (rho * betaS));
        end
        
        function val = getDensity(T,P)
           val = NitrousFluidCoolProp.getProperty(FluidProperty.DENSITY,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P); 
        end
        
        %Function to get the density of the non-saturated liquid (Kg/m^3)
        %at a given Temp (K) and Pressure (Pa)
        function val = getLiquidDensity(T,P)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.DENSITY,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        %Function to get the density of the non-saturated gas (Kg/m^3)
        %at a given Temp (K) and Pressure (Pa)
        function val = getGasDensity(T,P)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.DENSITY,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        function val = getSpecificEnthalpy(T,P)
            val = NitrousFluidCoolProp.getProperty(FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        function val = getGasSpecificEnthalpy(T,P)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        function val = getLiquidSpecificSaturationEnthalpy(T)
            %PVap = SaturatedNitrous.getVapourPressure(T);
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,T,FluidProperty.VAPOR_QUALITY,0);
        end
        
        function val = getEnthalpyOfVaporization(T)
            PVap = SaturatedNitrous.getVapourPressure(T);
            val = NitrousFluid.getGasSpecificEnthalpy(T,PVap) - NitrousFluid.getLiquidSpecificSaturationEnthalpy(T);
        end
        
        function val = getGasSpecificSaturationEnthalpy(T)
            %PVap = SaturatedNitrous.getVapourPressure(T);
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,T,FluidProperty.VAPOR_QUALITY,1);
        end

        %Function to get the specific enthalpy of the liquid in J/Kg
        function val = getLiquidSpecificEnthalpy(T,P)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        function Z = getCompressibilityFactor(T,P)
           Z = NitrousFluidCoolProp.getProperty(FluidProperty.COMPRESSIBILITY_FACTOR,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P); 
        end
        
        %Function to get the compressibility factor, Z, of the gas at a
        %given temperature and pressure. Is calculated by looking up known
        %density for this T and P and using Z=P/(rho)RT
        function Z = getGasCompressibilityFactor(T,P)
            Z = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.COMPRESSIBILIY_FACTOR,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
%             rho = NitrousFluid.getGasDensity(T,P);
%             Z = P / (rho * NitrousFluid.getGasConstant() * T);
        end
        
        %Function to get the compressibility factor, Z, of the liquid at a
        %given temperature and pressure. Is calculated by looking up known
        %density for this T and P and using Z=P/(rho)RT. Yes am using a
        %weird equation of state for this liquid (Making it match the gas for ease), but doesn't actually
        %really matter as this Z value is picked from real data to make it
        %match
        function Z = getLiquidCompressibilityFactor(T,P)
            Z = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.COMPRESSIBILITY_FACTOR,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
%             rho = NitrousFluid.getLiquidDensity(T,P);
%             Z = P / (rho * NitrousFluid.getGasConstant() * T);
        end
        
        function val = getCp(T,P)
            val = NitrousFluidCoolProp.getProperty(FluidProperty.CP,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);;
        end
        
        %Function to get the specific heat capacity at constant pressure for the
        %liquid (In J/K/Kg)
        %at a given Temp (K) and Pressure (Pa)
        function val = getLiquidCp(T,P)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.CP,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        %Function to get the isobaric coefficient of expansion (1/K)
        %at a given Temp (K) and Pressure (Pa)
        function val = getGasIsobaricCoeffOfExpansion(T,P)
            P1 = P/1000; %Need gas in kPa using tabulated data
            data = NitrousFluid.getDataFromFile(['nitrousRawData',filesep,'gasIsobaricCoefficientOfExpansion.txt']); 
            val = NitrousFluid.twoColInterp(data,1,2,3,T,P1);
        end
        
        %Function to get the isobaric coefficient of expansion (1/K)
        %at a given Temp (K) and Pressure (Pa)
        function val = getLiquidIsobaricCoeffOfExpansion(T,P)
            P1 = P/1000; %Need gas in kPa using tabulated data
            data = NitrousFluid.getDataFromFile(['nitrousRawData',filesep,'liquidIsobaricCoefficientOfExpansion.txt']); 
            val = NitrousFluid.twoColInterp(data,1,2,3,T,P1);
        end
        
        function val = getTemperatureForPressureEntropy(P,s)
            val = NitrousFluidCoolProp.getProperty(FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P,FluidProperty.SPECIFIC_ENTROPY,s);
        end
        
        %Function to get the temperature (K) of the gas for a given
        %pressure (Pa) and specific entropy (J/K/Kg)
        function val = getGasTemperatureForPressureEntropy(P,s)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P,FluidProperty.SPECIFIC_ENTROPY,s);
        end
        
        %Function to get the temperature (K) of the liquid for a given
        %pressure (Pa) and specific entropy (J/K/Kg)
        function val = getLiquidTemperatureForPressureEntropy(P,s)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P,FluidProperty.SPECIFIC_ENTROPY,s);
        end
        
        function val = getPressureForTemperatureEntropy(T,s)
           val = NitrousFluidCoolProp.getProperty(FluidProperty.PRESSURE,FluidProperty.TEMPERATURE,T,FluidProperty.SPECIFIC_ENTROPY,s);; 
        end
        
        %Function to get the pressure (Pa) of the gas for a given
        %temperature (T) and specific entropy (J/K/Kg)
        function val = getGasPressureForTemperatureEntropy(T,s)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.PRESSURE,FluidProperty.TEMPERATURE,T,FluidProperty.SPECIFIC_ENTROPY,s);
        end
        
        %Function to get the pressure (Pa) of the liquid for a given
        %temperature (T) and specific entropy (J/K/Kg)
        function val = getLiquidPressureForTemperatureEntropy(T,s)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.PRESSURE,FluidProperty.TEMPERATURE,T,FluidProperty.SPECIFIC_ENTROPY,s);
        end
        
        function val = getPressureForTemperatureEnthalpy(T,h)
            val = NitrousFluidCoolProp.getProperty(FluidProperty.PRESSURE,FluidProperty.TEMPERATURE,T,FluidProperty.SPECIFIC_ENTHALPY,h);
        end

        %Function to get the pressure (Pa) of the gas for a given specific
        %enthalpy (J/Kg) and temperature (K)
        function val = getGasPressureForTemperatureEnthalpy(T,h)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.PRESSURE,FluidProperty.TEMPERATURE,T,FluidProperty.SPECIFIC_ENTHALPY,h);
        end
        
        function val = getTemperatureForPressureEnthalpy(P,h)
            try
                val = NitrousFluidCoolProp.getProperty(FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P,FluidProperty.SPECIFIC_ENTHALPY,h);
            catch
                val = NitrousFluidCoolProp.getProperty(FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P,FluidProperty.SPECIFIC_ENTHALPY,h+2);
            end
        end
        
        %Function to get the temperature (K) of the gas for a given specific
        %enthalpy (J/Kg) and pressure (Pa)
        function val = getGasTemperatureForPressureEnthalpy(P,h)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P,FluidProperty.SPECIFIC_ENTHALPY,h);
        end
        
        %Function to get the pressure (Pa) of the liquid for a given specific
        %enthalpy (J/Kg) and temperature (K)
        function val = getLiquidPressureForTemperatureEnthalpy(T,h)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.PRESSURE,FluidProperty.TEMPERATURE,T,FluidProperty.SPECIFIC_ENTHALPY,h);
        end
        
        %Function to get the temperature (K) of the gas for a given specific
        %enthalpy (J/Kg) and pressure (Pa)
        function val = getLiquidTemperatureForPressureEnthalpy(P,h)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.TEMPERATURE,FluidProperty.PRESSURE,P,FluidProperty.SPECIFIC_ENTHALPY,h);
        end
        
        function Cv = getCv(T,P)
            Cv = NitrousFluidCoolProp.getProperty(FluidProperty.CV,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        %Function to get the specific heat capacity at constant volume for the
        %liquid (In J/K/Kg)
        %at a given Temp (K) and Pressure (Pa)
        function Cv = getLiquidCv(T,P)
            Cv = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.CV,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        %Function to get the specific heat capacity at constant pressure for the
        %gas (In J/K/Kg)
        %at a given Temp (K) and Pressure (Pa)
        function val = getGasCp(T,P)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.CP,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        %Function to get the specific heat capacity for the saturated mixture (In J/K/Kg)
        %at a given Temp (K) and Pressure (Pa)
        function val = getSaturatedHeatCapacity(T)
            data = NitrousFluid.getDataFromFile(['nitrousRawData',filesep,'saturationHeatCapacity.txt'],3); 
            val = NitrousFluid.oneColInterp(data,1,2,T);
            %Val is in J/K/mol, convert to J/K/Kg
            val = val / NitrousFluid.getMolarMass();
        end
        
        function val = getPressureFromDensityEnthalpy(rho,h)
            data = NitrousFluid.getDataFromFile(['preBakedData',filesep,'pressureDensityEnthalpy.mat']);
            val = NitrousFluid.twoColInterp(data,2,1,3,h,rho);
        end
        
        function val = getLiquidEntropy(T,P)
             val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_ENTROPY,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        %Function to get the entropy of the liquid in equilibrium with the
        %gas in a saturated state
        function val = getSaturatedLiquidEntropy(T)
            P = SaturatedNitrous.getVapourPressure(T);
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.SPECIFIC_ENTROPY,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        function val = getPressure(T,density)
           val = NitrousFluidCoolProp.getProperty(FluidProperty.PRESSURE,FluidProperty.TEMPERATURE,T,FluidProperty.DENSITY,density);
        end

        %Function to get the pressure (Pa) of the liquid for a given T (K)
        %and density (kg/m^3)
        function val = getLiquidPressure(T,density)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.LIQUID,FluidProperty.PRESSURE,FluidProperty.TEMPERATURE,T,FluidProperty.DENSITY,density);
        end
        
        %Function to get the pressure (Pa) of the gas for a given T (K)
        %and density (kg/m^3)
        function val = getGasPressure(T,density)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.PRESSURE,FluidProperty.TEMPERATURE,T,FluidProperty.DENSITY,density);
        end
        
        function val = getEntropy(T,P)
           val = NitrousFluidCoolProp.getProperty(FluidProperty.SPECIFIC_ENTROPY,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P); 
        end
        
        function val = getGasEntropy(T,P)
            val = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.SPECIFIC_ENTROPY,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        %Function to get the specific heat capacity at constant volume for the
        %gas (In J/K/Kg)
        %at a given Temp (K) and Pressure (Pa)
        function Cv = getGasCv(T,P)
            Cv = NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.CV,FluidProperty.TEMPERATURE,T,FluidProperty.PRESSURE,P);
        end
        
        function val = getSpecificHeatRatio(T,P)
           val = NitrousFluid.getCp(T,P) / NitrousFluid.getCv(T,P); 
        end
        
        %Function to get the specific heat ratio (Cp/Cv) for the liquid at
        %a given Temp (K) and Pressure (Pa)
        function val = getLiquidSpecificHeatRatio(T,P)
            val = NitrousFluid.getLiquidCp(T,P) / NitrousFluid.getLiquidCv(T,P);
        end
        
        %Function to get the specific heat ratio (Cp/Cv) for the gas at
        %a given Temp (K) and Pressure (Pa)
        function val = getGasSpecificHeatRatio(T,P)
            val = NitrousFluid.getGasCp(T,P) / NitrousFluid.getGasCv(T,P);
        end
    end
end