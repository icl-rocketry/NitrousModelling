clear
clc
close all
%Diameter of the bore (hole cut through ball) in m
boreDiameter = 4e-3; %Arbitrary value
%Diameter of the ball within the valve in m
ballDiameter = 6e-3; %Arbitrary value
%openAmt = 0->1 for closed=0 and open=1
getArea = @(openAmt) FullBoreBallValve.getOrificeArea(openAmt.*(pi/2),boreDiameter,ballDiameter);
openAmts = linspace(0,1,20);
areas = getArea(openAmts);
plot(openAmts, areas.*1000000);
title("Full bore ball valve area seen by flow, bore="+boreDiameter+", ball="+ballDiameter);
xlabel("Opening angle fraction");
ylabel("Area (mm^2)");