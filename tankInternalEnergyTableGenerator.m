clear;
% T = 182.23:1:SaturatedNitrous.T_CRIT;
EMax = 4.2e5; %SPECIFIC internal energy (eg per unit mass of nitrous)
E = 0:500:EMax;
u = 0;
mNitrous = 0;
mNitrousMax = 20; %Kg
mStep = 0.010;
filename='preBakedData/nitrousTankInternalEnergiesv3.mat';

disp("Creating data struts...");
drawnow;
tankCopy = GeometricNitrousTank(20+273.15,mNitrous,1,0.01767,0.8);
tankVol = tankCopy.tankTotalVolume;
disp(tankVol);
pause(2);
for i=1:length(E)
    mNitrous = 0;
    disp(((E(i)*100)/EMax)+"%");
    drawnow;
    while (mNitrous <= mNitrousMax)
        %GeometricNitrousTank(temp,mNitrousTotal,tankHeight,crossSectionalArea,gasVentHoleHeight)
%         if(tankCopy.mLiquid < 0 || tankCopy.mVapour < 0)
%             %If not saturated, don't tabulate it as value is probably wrong
%             %disp("Skipping "+T(i)+" "+mNitrous);
%             mNitrous = mNitrous + mStep;
%             continue;
%         end
        u = u+1;
        %Calculate enthalpy for this T and P
%         h = FluidType.NITROUS_LIQUID.getSpecificEnthalpy(T(i),P);
%         str = sprintf('%d \t %d \t %d',T(i),P,h);
%         disp(str);
%         drawnow;
        data{1,u} = E(i);
        data{2,u} = mNitrous / tankVol;
        disp(mNitrous+" -> "+data{2,u});
        data{3,u} = {};
        
        %Advance m nitrous
        mNitrous = mNitrous + mStep;
    end
end

disp("Created data struts ("+length(data)+" data points), STARTING calculations...");
drawnow;

if(true)
    return;
end

dataMap = containers.Map('KeyType','char','ValueType','any');

dataCopy = data;
% lines = cell(length(data));
% fileID = fopen('preBakedData/saturatedTankInternalEnergy.txt','w');
% fprintf(fileID,['Table: Internal Energy [Saturated Tank]\r\n']);
% fprintf(fileID,['Temperature {K} \t Total Nitrous Mass Per Nitrous Tank Vol {Kg/m^3} \t Specific Internal Energy (Per Kg Total Nitrous) {J/Kg}\r\n']);
len = length(data);
tankVol = tankCopy.tankTotalVolume;
tic;
parfor k=1:len
    drawnow;
    if ~isempty(data{3,k})
       continue; 
    end
    mNitrous = dataCopy{2,k} * tankVol;
    E = dataCopy{1,k} * mNitrous;
    tank = GeometricNitrousTank(20+273.15,mNitrous,1,0.01767,0.8);
    try
        tank.setNitrousMassEnergyAndAdjustTemp(mNitrous,E);
    catch excep
        %disp(excep);
        data{3,k} = -1;
        continue; %Invalid datapoint
    end
    data{3,k} = tank.temp;
end
toc;
disp("Finished calculating internal energies correlation to temp, writing to map then saving...");
drawnow;
for z=1:length(data)
    specificE = dataCopy{1,z};
    nitrousPerVol = dataCopy{2,z};
    key = [num2str(specificE),'|',num2str(nitrousPerVol)];
    if(data{3,z} ~= -1)
        dataMap(key) = data{3,z};
    end
end
disp("Saving map....");
drawnow;

save(filename,'dataMap');
disp("Done!");
drawnow;
