classdef SystemVisualization < matlab.mixin.Copyable
    properties(Access = private)
        figure; 
        axis;
        texts;
    end
    methods
        function obj = SystemVisualization()
           obj.figure = figure(); 
           axis([0 125 0 100]); %Creates an axis
           obj.axis = obj.figure.CurrentAxes;
           obj.drawFixed();
        end
        
        function drawText(obj,i,x,y,txt)
           if i<=length(obj.texts)
               try 
                    delete(obj.texts(i)); 
               catch
                   
               end
           end
           obj.texts(i) = text(obj.axis,x,y,txt);
        end
        
        function drawFixed(obj)
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
        end
        
        function draw(obj, x, u,t)
            %Draw stuff that changes
            obj.drawText(1,46,45,['\theta_1 = ',num2str(u(1))]);
            %text(obj.axis,50,45,['\theta_1 = ',num2str(u(1))]);
            obj.drawText(2,1,70,['\theta_2 = ',num2str(u(2))]);
            obj.drawText(3,50,90,['T_a = ',num2str(x(1)-273.15),' C']);
            obj.drawText(4,76,60,['T_e = ',num2str(x(2)-273.15),' C']);
            obj.drawText(5,76,50,['m_e = ',num2str(x(3)),' Kg']);
            obj.drawText(6,21,60,['T_i = ',num2str(x(4)-273.15),' C']);
            obj.drawText(7,21,50,['m_i = ',num2str(x(5)),' Kg']);
            obj.drawText(8,101,70,['Q = ',num2str(u(3)),' W']);
            obj.drawText(9,101,40,['T_f = ',num2str(x(6)-273.15),' C']);
            
            title(obj.axis,['t = ',num2str(t),' s']);
        end
    end
end