classdef FullBoreBallValve
    methods (Static,Access=public)
        function orifice = getOrificeArea(theta,boreDiameter,ballDiameter)
            o = @(z) FullBoreBallValve.getOpeningLength(theta,boreDiameter,ballDiameter,z);
            orifice = integral(o,0,ballDiameter,'ArrayValued',true);
        end

        %Get the length from the tip of the bore to the housing
        %for a FULL BORE ball valve at opening angle theta in radians (theta=pi/2 is fully open)
        %with a given bore diameter and ball diameter at a given height, z, from
        %the base of the ball. Eg. integrate this wrt z from 0 to <ball diameter> to get the orifice area
        function opening = getOpeningLength(theta,boreDiameter,ballDiameter,z)
            z = z-(ballDiameter/2); %Shift z coordinate to match what used for derivation

            %Radius of circlular section of ball valve (circle looking top down) at this height z
            R = 0.5.*sqrt((ballDiameter.^2) - 4.*(z.^2));
            %Bore diameter at this circular section of ball valve (circle looking
            %down) at this height z - eg. the length of the chord on that circle
            c = sqrt((boreDiameter.^2) - 4.*(z.^2)) + (theta - theta); %+ and - theta is a hack to get array dimens to match for array input

            %Arg is the angle formed between the centerline of the bore and the
            %casing - a bit of a fun derivation to do. Look on the drive to see how
            %got here.
            arg = (theta/2) - (pi/4) + 0.5*acos(1-((c.^2)./(2.*(R.^2))));
            %The opening is 2rsin(arg)
            opening = 2.*R.*sin(arg);
            for i=1:length(opening)
               if(opening(i) <= 0 || isnan(opening(i)) || imag(c(i)) ~= 0)
                  opening(i) = 0; 
               end
            end
        end
    end
end

