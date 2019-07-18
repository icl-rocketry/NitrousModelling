clear
clc

%Ball Valve (Data for generic full bore) Cv vs angle function finder
amtOpen = [0 10 20 30 40 50 60 70 80 90] ./ 90;
typicalValveCvFrac = [0 3 5 9 15 23 35 58 77 90] ./ 90;
plot(amtOpen,typicalValveCvFrac);
P = polyfit(amtOpen,typicalValveCvFrac,6);
hold on;
fplot(@(x)polyval(P,x),[0,1]);
legend('Raw data', 'Fitted');
zeroPt = polyval(P,0);
disp(zeroPt);
P(length(P)) = P(length(P)) - zeroPt;
disp(polyval(P,0));
