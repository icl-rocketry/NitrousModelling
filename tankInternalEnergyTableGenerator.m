clear;
T = 240:1:SaturatedNitrous.T_CRIT;
u = 0;
mNitrous = 0;
mNitrousMax = 60; %Kg
mStep = 0.1;
for i=1:length(T)
    mNitrous = 0;
    while (mNitrous <= mNitrousMax)
        tankCopy = GeometricNitrousAndAirTank(T(i),mNitrous,0,1,0.01767,0.8);
        if(tankCopy.mLiquid < 0 || tankCopy.mVapour < 0)
            %If not saturated, don't tabulate it as value is probably wrong
            %disp("Skipping "+T(i)+" "+mNitrous);
            mNitrous = mNitrous + mStep;
            continue;
        end
        u = u+1;
        %Calculate enthalpy for this T and P
%         h = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(T(i),P);
%         str = sprintf('%d \t %d \t %d',T(i),P,h);
%         disp(str);
%         drawnow;
        data(u,1) = T(i);
        data(u,2) = mNitrous / tankCopy.volNitrous;
        data(u,3) = -1;
        
        %Advance m nitrous
        mNitrous = mNitrous + mStep;
    end
end

dataCopy = data;
lines = cell(length(data));
fileID = fopen('saturatedTankInternalEnergy.txt','w');
fprintf(fileID,['Table: Internal Energy [Saturated Tank]\r\n']);
fprintf(fileID,['Temperature {K} \t Total Nitrous Mass Per Nitrous Tank Vol {Kg/m^3} \t Specific Internal Energy (Per Kg Total Nitrous) {J/Kg}\r\n']);
parfor k=1:length(data)
    T = dataCopy(k,1);
    mNitrous = dataCopy(k,2) * tankCopy.volNitrous;
    tank = GeometricNitrousAndAirTank(T,mNitrous,0,1,0.01767,0.8);
    E = tank.getInternalEnergy() / mNitrous;
    data(k,3) = E;
    lines(k) = {sprintf('%d \t %d \t %d',T,(mNitrous/tankCopy.volNitrous),E)};
end
disp("Finished calculating internal energies");
drawnow;

for i=1:length(lines)
   fprintf(fileID,[cell2mat(lines(i)),'\r\n']); 
end
fclose(fileID);
