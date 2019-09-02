clear
clc
close all

vis = SystemVisualization();
x = [273.15 + 25
    273.15 + 30;
    10;
    273.15 + 20;
    2];
u = [0;
    0;
    0];
t = 0;
vis.draw(x,u,t);
while (t < 100)
   pause(0.01);
   t = t+0.01;
   vis.draw(x,u,t);
end