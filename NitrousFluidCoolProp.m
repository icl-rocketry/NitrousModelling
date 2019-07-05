classdef NitrousFluidCoolProp
   properties(Constant)
        FLUID_NAME = 'NitrousOxide';
    end
    
   methods (Static,Access=public)
       %Get a property of the fluid in the given phase at a given state defined by the args
       %passed.
       %Eg. NitrousFluidCoolProp.getPropertyForPhase(FluidPhase.GAS,FluidProperty.DENSITY,FluidProperty.PRESSURE,50e5,FluidProperty.TEMPERATURE,300)
       function val = getPropertyForPhase(fluidPhase,propertyToGet,varargin)
           if(isempty(varargin))
              val = NitrousFluidCoolProp.getProperty(propertyToGet);
              return;
           end
           
           args2 = varargin;
           firstArg = args2{1}(1);
           if isa(firstArg,'FluidProperty')
               firstArg = firstArg.getCoolPropInput();
           end
           if ~isa(firstArg,'char')
               error('First arg defining fluid state must be either a FluidProperty or a char sequence');
           end
           firstArg = strcat([firstArg,'|',fluidPhase.getCoolPropInput()]);
           args2{1} = firstArg;
           val = NitrousFluidCoolProp.getProperty(propertyToGet,args2);
       end
       
        %Get a property of the fluid at a given state defined by args
        %passed.
        %Eg. NitrousFluidCoolProp.getProperty(FluidProperty.DENSITY,FluidProperty.PRESSURE,50e5,FluidProperty.TEMPERATURE,300)
        function val = getProperty(propertyToGet,varargin)
            if isa(propertyToGet,'FluidProperty')
                propToGetInput = propertyToGet.getCoolPropInput();
            else 
                propToGetInput = propertyToGet;
            end
            cmd = strcat(['py.CoolProp.CoolProp.PropsSI(''',propToGetInput,''',']);
            if length(varargin)==1 && isa(varargin,'cell')
                arg1 = varargin{1};
                if isa(arg1,'cell')
                   varargin = arg1;
                end
            end
            if (rem(length(varargin),2)~=0)
               error('Properties specified to define the fluid state not fully defined'); 
            end
            for i=1:length(varargin)
               arg = varargin{i};
               if isa(arg,'FluidProperty')
                   arg = arg.getCoolPropInput();
               end
               if isa(arg,'char')
                   cmd = strcat([cmd,'''',arg,''',']);
               else
                   cmd = strcat([cmd,num2str(arg),',']);
               end
            end
            
            cmd = strcat([cmd,'''',NitrousFluidCoolProp.FLUID_NAME,''')']);
            val = eval(cmd);
%             val = py.CoolProp.CoolProp.PropsSI(propertyToGet.getCoolPropInput(),'NitrousOxide');
        end
   end
end