clear;
T = 182.33:1:309.521;
% P = 55e5;
PMax = 10000e3;
PStep = 200e3;
j = 0;
for i=1:length(T)
    PMin = SaturatedNitrous.getVapourPressure(T(i));
    P = PMin;
    while (P <= PMax)
        j = j+1;
        %Calculate enthalpy for this T and P
%         h = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(T(i),P);
%         str = sprintf('%d \t %d \t %d',T(i),P,h);
%         disp(str);
%         drawnow;
        data(j,1) = T(i);
        data(j,2) = P;
        data(j,3) = -1;
        
        %Advance P
        if(P == PMin)
            P = ceil(PMin / PStep) * PStep;
            if(P == PMin)
               P = P + PStep; 
            end
        else
            P = P + PStep;
        end
    end
end
dataCopy = data;
lines = cell(length(data));
fileID = fopen('nitrousLiquidInternalEnergy.txt','w');
fprintf(fileID,['Table: Internal Energy [L]\r\n']);
fprintf(fileID,['Temperature {K} \t Pressure {Pa} \t Specific Internal Energy {J/Kg}\r\n']);
parfor k=1:length(data)
    T = dataCopy(k,1);
    P = dataCopy(k,2);
    h = FluidType.NITROUS_LIQUID.getSpecificInternalEnergy(T,P);
    data(k,3) = h;
    lines(k) = {sprintf('%d \t %d \t %d',T,P,h)};
end
disp("Finished calculating internal energies");
drawnow;

for i=1:length(lines)
   fprintf(fileID,[cell2mat(lines(i)),'\r\n']); 
end
fclose(fileID);
