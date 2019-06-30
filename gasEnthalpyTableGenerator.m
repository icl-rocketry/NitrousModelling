clear;
T = 197:1:400;
% P = 55e5;
PMax = 10000e3;
PStep = 200e3;
j = 0;
for i=1:length(T)
    VP = SaturatedNitrous.getVapourPressure(T(i));
    if(T(i) > SaturatedNitrous.T_CRIT)
       VP = 10000e3; %Ignore maximum after critical pt 
    end
    PMax = min(10000e3,VP);
    P = 200e3;
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
        P = P + PStep;
    end
end
dataCopy = data;
lines = cell(length(data));
fileID = fopen('nitrousGasEnthalpy.txt','w');
fprintf(fileID,['Table: Enthalpy [G]\r\n']);
fprintf(fileID,['Temperature {K} \t Pressure {Pa} \t Specific Enthalpy {J/Kg}\r\n']);
parfor k=1:length(data)
    T = dataCopy(k,1);
    P = dataCopy(k,2);
    h = FluidType.NITROUS_GAS.getSpecificEnthalpy(T,P);
    data(k,3) = h;
    lines(k) = {sprintf('%d \t %d \t %d',T,P,h)};
end
disp("Finished calculating enthalpies");
drawnow;

for i=1:length(lines)
   fprintf(fileID,[cell2mat(lines(i)),'\r\n']); 
end
fclose(fileID);
