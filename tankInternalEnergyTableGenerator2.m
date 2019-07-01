clear;
mNitrousRange = 0:0.1:60; %Kg
u = 0;
TMax = SaturatedNitrous.T_CRIT;
TStep = 1;
for i=1:length(mNitrousRange)
    T = 240;
    while (T <= TMax)
        tankCopy = GeometricNitrousAndAirTank(T,mNitrousRange(i),0,1,0.01767,0.8);
        if(tankCopy.mLiquid < 0 || tankCopy.mVapour < 0)
            %If not saturated, don't tabulate it as value is probably wrong
            %disp("Skipping "+T(i)+" "+mNitrous);
            T = T + TStep;
            continue;
        end
        u = u+1;
        %Calculate enthalpy for this T and P
%         h = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(T(i),P);
%         str = sprintf('%d \t %d \t %d',T(i),P,h);
%         disp(str);
%         drawnow;
        data(u,1) = mNitrousRange(i) / tankCopy.volNitrous;
        data(u,2) = T;
        data(u,3) = -1;
        
        %Advance T
        T = T + TStep;
    end
end

dataCopy = data;
lines = cell(length(data));
fileID = fopen('saturatedTankInternalEnergy2.txt','w');
fprintf(fileID,['Table: Internal Energy [Saturated Tank]\r\n']);
fprintf(fileID,['Total Nitrous Mass Per Nitrous Tank Vol {Kg/m^3} \t Temperature {K} \t Specific Internal Energy (Per Kg Total Nitrous) {J/Kg}\r\n']);
parfor k=1:length(data)
    T = dataCopy(k,2);
    mNitrous = dataCopy(k,1) * tankCopy.volNitrous;
    tank = GeometricNitrousAndAirTank(T,mNitrous,0,1,0.01767,0.8);
    E = tank.getInternalEnergy() / mNitrous;
    data(k,3) = E;
    lines(k) = {sprintf('%d \t %d \t %d',(mNitrous/tankCopy.volNitrous),T,E)};
end
disp("Finished calculating internal energies");
drawnow;

for i=1:length(lines)
   fprintf(fileID,[cell2mat(lines(i)),'\r\n']); 
end
fclose(fileID);
