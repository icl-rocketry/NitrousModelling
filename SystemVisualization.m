classdef SystemVisualization < matlab.mixin.Copyable
    properties
        figure; 
        axis;
    end
    methods
        function obj = SystemVisualization()
           obj.figure = figure(); 
           axis([0 125 0 100]); %Creates an axis
           obj.axis = obj.figure.CurrentAxes;
        end
        
        function draw(obj, x, u,t)
            cla(obj.axis); %Clear axis
            %Draw left (internal) tank at (20,20) width 25 and height 60
            rectangle(obj.axis,'Position',[20 20 25 60]);
            %Draw right (external) tank at (75,20) width 25 and height 60
            rectangle(obj.axis,'Position',[75 20 25 60]);
            
            %Draw pipe between tanks
            rectangle(obj.axis,'Position',[45 40 30 10]);
            %Draw vent pipe
            rectangle(obj.axis,'Position',[0 65 20 10]);
            
            %Labels
            text(obj.axis,25,75,'Int. Tank');
            text(obj.axis,80,75,'Ext. Tank');
            text(obj.axis,50,45,['\theta_1 = ',num2str(u(1))]);
            text(obj.axis,5,70,['\theta_2 = ',num2str(u(2))]);
            text(obj.axis,50,90,['T_a = ',num2str(x(1)-273.15),' C']);
            text(obj.axis,76,60,['T_e = ',num2str(x(2)-273.15),' C']);
            text(obj.axis,76,50,['m_e = ',num2str(x(3)),' Kg']);
            text(obj.axis,21,60,['T_i = ',num2str(x(4)-273.15),' C']);
            text(obj.axis,21,50,['m_i = ',num2str(x(5)),' Kg']);
            text(obj.axis,80,85,['Q = ',num2str(u(3)),' W']);
            
            title(obj.axis,['t = ',num2str(t),' s']);
        end
    end
end