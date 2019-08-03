clear;
%P = 200e3:200e3:100e5;
rho = 1:2:1200;
% P = 55e5;
enthalpy = -2.5e5:0.01e5:5.2e5;
j = 0;
tic; P2 = findP(1000,1e5); toc;
for i=1:length(enthalpy)
    for k=1:length(rho)
        j = j+1;
        data(j,1) = rho(k);
        data(j,2) = enthalpy(i);
        data(j,3) = -1;
    end
end
dataCopy = data;
lines = cell(size(data));
fileID = fopen('pressureDensityEnthalpy.txt','w');
fprintf(fileID,['Table: Enthalpy [All Phase]\r\n']);
fprintf(fileID,['Density {Kg/m^3} \t Specific Enthalpy {J/Kg} \t Pressure {Pa}\r\n']);
parfor k=1:length(data)
    rho = dataCopy(k,1);
    h = dataCopy(k,2);
    %hErr = @(P) h - NitrousFluidCoolProp.getProperty(FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.PRESSURE,P,FluidProperty.DENSITY,rho);
    %P2 = fzero(hErr,30e5);
    try
        P2 = findP(rho,h);
    catch
        P2 = -1;
    end
    %h = FluidType.NITROUS_GAS.getSpecificEnthalpy(T,P);
    data(k,3) = P2;
    lines(k) = {sprintf('%d \t %d \t %d',rho,h,P2)};
end
disp("Finished calculating enthalpies, writing to file");
drawnow;

for i=1:length(lines)
   fprintf(fileID,[cell2mat(lines(i)),'\r\n']); 
end
fclose(fileID);

disp("Finished writing to file");
drawnow;

function P2 = findP(rho,h)
    P2 = abs(real(fzero(@(P) hErrFunc(P),30e5)));
    function hErr = hErrFunc(P)
        P = abs(real(P));
        hErr = h - NitrousFluidCoolProp.getProperty(FluidProperty.SPECIFIC_ENTHALPY,FluidProperty.PRESSURE,P,FluidProperty.DENSITY,rho);
    end
end
