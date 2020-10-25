clear
clc

launchTempC = 28; %Corresponds to 60.38 bar pressure
ullage = 0.2; %20% ullage
cylinderTankDiam = 150e-3; %150mm internal diameter
fractionNitrousRemainingExtra = 0.1; %Have 10% of required nitrous sat at tank at end of flight
requiredOxidiserMassForEngine = 6; %Kg

launchTemp = launchTempC + 273.15;
tankCrossSection = 0.25*pi*cylinderTankDiam^2;
hReq = sizeMeATank(launchTemp,requiredOxidiserMassForEngine,tankCrossSection,ullage,fractionNitrousRemainingExtra);
vReq = hReq * tankCrossSection;

disp("Height of simple cylinder tank: "+hReq+" m");
disp("Internal volume of simple cylinder tank: "+vReq+" m^3");
