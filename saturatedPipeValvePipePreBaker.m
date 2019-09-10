clear
clc
%For gas on saturation line
filename = 'preBakedData/saturatedGasPipeValveFlowRatesSmallerValve.mat';
%Gas one calculated with internal diameter of 4.8mm
%Liquid one calculated with internal diamater of 10.2108mm
pipeInternalDiameter = 4.8e-3;
valveFullyOpenFlowCoefficient = 0.2; %12 for 1/4in and 0.2 for 1/8in
upstreamQuality = 1; %1 is vapour, 0 is liquid

pipe1 = FluidPipe(0.25*pi*(pipeInternalDiameter).^2,1);
valveOpenAmt = 0:0.025:1;
% valve = BallValve(12,valveOpenAmt(1));
% disp("Flow coeff: "+valve.getFlowCoefficient());
pipe2 = FluidPipe(0.25*pi*(pipeInternalDiameter).^2,1);
pressures = 1e5:100e3:72e5;

disp("Generating data structure...");
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
load('preBakedData/dataNotMapPartial.mat','dataNotMap');
dataCopy = dataNotMap;
data = containers.Map('KeyType','char','ValueType','any');

disp("Starting calculations...");
drawnow;
len = length(dataNotMap);
tic;
parfor z=1:len
    upstreamPressure = dataCopy{1,z};
    downstreamPressure = dataCopy{2,z};
    if ~isempty(dataNotMap{3,z})
       continue; 
    end
    upstreamTemp = SaturatedNitrous.getSaturationTemperature(upstreamPressure);
    mdot = zeros(1,length(valveOpenAmt));
    for i=1:length(valveOpenAmt)
        valve = BallValve(valveFullyOpenFlowCoefficient,valveOpenAmt(i));
        pvp = PipeValvePipe(pipe1,valve,pipe2);
        try
            [~,mdot(i),~,~] = pvp.getDownstreamTemperatureMassFlowFromPressureChange(downstreamPressure-upstreamPressure,FluidType.NITROUS_GENERAL,upstreamTemp,upstreamPressure,upstreamQuality,0);
        catch excep
            disp("Upstream P: "+upstreamPressure);
            disp("Downstream P: "+downstreamPressure);
            disp("Valve open amt: "+valveOpenAmt(i));
            drawnow;
            rethrow(excep);
        end
    end
    mdotPolynomialFitCoeffs = polyfit(valveOpenAmt,mdot,10);
    dataNotMap{3,z} = mdotPolynomialFitCoeffs;    
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

save(filename,'data');
disp("Done!");
drawnow;