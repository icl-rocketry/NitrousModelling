clear;
P = 200e3:200e3:10000e3;
% P = 55e5;
TAbsMin = 150;
TAbsMax = 340;
TStep = 1;
j = 0;
for i=1:length(P)
    TMin = SaturatedNitrous.getSaturationTemperature(P(i));
    if(TMin == -1) %Fluid super critical, liquid not exist
       TMin = SaturatedNitrous.T_CRIT;
    end
    T = TMin;
    while (T <= TAbsMax)
        j = j+1;
        %Calculate enthalpy for this T and P
%         h = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(T(i),P);
%         str = sprintf('%d \t %d \t %d',T(i),P,h);
%         disp(str);
%         drawnow;
        data(j,1) = P(i);
        data(j,2) = T;
        data(j,3) = -1;
        
        %Advance T
        if(T == TMin)
            T = ceil(TMin / TStep) * TStep;
            if(T == TMin)
               T = T + TStep; 
            end
        else
            T = T + TStep;
        end
    end
end
dataCopy = data;
lines = cell(length(data));
fileID = fopen('nitrousGasInternalEnergy2.txt','w');
fprintf(fileID,['Table: Internal Energy [G]\r\n']);
fprintf(fileID,['Pressure {Pa} \t Temperature {K} \t Specific Internal Energy {J/Kg}\r\n']);
parfor k=1:length(data)
    T = dataCopy(k,2);
    P = dataCopy(k,1);
    h = FluidType.NITROUS_GAS.getSpecificInternalEnergy(T,P);
    data(k,3) = h;
    lines(k) = {sprintf('%d \t %d \t %d',P,T,h)};
end
disp("Finished calculating internal energies");
drawnow;

for i=1:length(lines)
   fprintf(fileID,[cell2mat(lines(i)),'\r\n']); 
end
fclose(fileID);
