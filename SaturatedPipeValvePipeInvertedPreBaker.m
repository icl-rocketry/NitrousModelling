clear
clc
filenameIn = 'preBakedData/saturatedLiquidPipeValveFlowRates.mat';
filenameOut = 'preBakedData/saturatedLiquidPipeValveFlowRatesInverted.mat';

disp("Generating data structure...");
pressures = 1e5:100e3:72e5;
drawnow;
j=0;
for i=1:length(pressures)
    upstreamPressure = pressures(i);
    downstreamPressures = 1e5:100e3:upstreamPressure;
    for y=1:length(downstreamPressures)
        j = j+1;
        dataNotMap{1,j} = upstreamPressure;
        dataNotMap{2,j} = downstreamPressures(y);
        dataNotMap{3,j} = {};
    end
end
% load('preBakedData/dataNotMapPartial.mat','dataNotMap');
dataCopy = dataNotMap;
data = containers.Map('KeyType','char','ValueType','any');

disp("Starting calculations...");
drawnow;
len = length(dataNotMap);
tic;
parfor z=1:len
    spvp = SaturatedPipeValvePipeFast(filenameIn); %Matlab forces you to load this per thread...
    upstreamPressure = dataCopy{1,z};
    downstreamPressure = dataCopy{2,z};
    if ~isempty(dataNotMap{3,z})
       continue; 
    end
    valvePos = 0:0.0025:1;
    mdots = zeros(1,length(valvePos));
    for i=1:length(valvePos)
        mdots(i) = spvp.getMassFlow(upstreamPressure,downstreamPressure,valvePos(i));
    end
    valvePosPolynomialFitCoeffs = polyfit(mdots,valvePos,10);
    valvePosPolynomialFitCoeffs(length(valvePosPolynomialFitCoeffs)) = 0;
    dataNotMap{3,z} = valvePosPolynomialFitCoeffs;    
end
toc;
disp("Finished calculations, writing to map and then saving...");
drawnow;
for z=1:length(dataNotMap)
    upstreamPressure = dataCopy{1,z};
    downstreamPressure = dataCopy{2,z};
    key = [num2str(upstreamPressure),'|',num2str(downstreamPressure)];
    data(key) = dataNotMap{3,z}; 
end
disp("Saving map....");
drawnow;
% 
% key = [num2str(upstreamPressure),'|',num2str(downstreamPressure)];
%     data(key) = mdotPolynomialFitCoeffs;

save(filenameOut,'data');
disp("Done!");
drawnow;