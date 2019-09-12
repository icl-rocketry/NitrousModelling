classdef FillingSystemTrajectory
   properties
      dtStep;
      u; %Array of u at each value of t from t=0 onwards in steps of dt. Where each timestep's u is a vertical vector
      x; %Array of x at each value of t from t=0 onwards in steps of dt
   end
   methods(Static)
       function traj = concatTrajectoriesWithSameDt(traj1,traj2)
           dtStep = traj1.dtStep;
           x = [traj1.x,traj2.x];
           u = [traj1.u,traj2.u];
           traj = FillingSystemTrajectory(dtStep,x,u);
       end
   end
   methods
       function obj = FillingSystemTrajectory(dtStep,x,u)
           obj.dtStep = dtStep;
           obj.x = x;
           obj.u = u;
       end
       
       function [x,u] = getXAndUForTime(obj,t)
           if(t > obj.dtStep*(length(obj.x)-1))
               %t is greater than trajectory specifies so just return last
               %datapoint
               x = obj.x(:,length(obj.x));
               u = obj.u(:,length(obj.u));
               return;
           end
           
           i = floor(t/obj.dtStep); %Index of closest datapoint in discretised trajectory saved
           i2 = i+1; %Find two existing points either side
           i3 = i-1;
           if(i3 < 1)
              i3 = 1;
              i = 2;
              i2 = 3;
           end
           if(i2 > length(obj.x))
              i2 = length(obj.x);
              i = length(obj.x) - 1;
              i3 = length(obj.x) - 2;
              if(i3 < 1)
                 i3 = 1; 
              end
              if(i < 1)
                  i = 1;
              end
           end
           
           ti = (i-1)*obj.dtStep;
           ti2 = ti+obj.dtStep;
           ti3 = ti-obj.dtStep;
           
           %Get values and interpolate
           xi = obj.x(:,i);
           xi2 = obj.x(:,i2);
           xi3 = obj.x(:,i3);
           ui = obj.u(:,i);
           ui2 = obj.u(:,i2);
           ui3 = obj.u(:,i3);
           
           x = spline([ti3,ti,ti2],[xi3,xi,xi2],t);
           u = spline([ti3,ti,ti2],[ui3,ui,ui2],t);
       end
   end
    
end