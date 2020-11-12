%1D solver that uses newton raphson but is also bracketed
% - aims for quick convergence without newton raphson hitting a
% weird point and suddenly jumping miles away
% Root MUST lie between minX and maxX to work and f(minX) and f(maxX)
% should have a change of sign
% Created whilst trying reliably solve a function in
% GeometricNitrousAndAirTank since fzero wasn't able to solve it
%
% By Eddie Brown
function [x] = betterfzero(func,x0,minX,maxX,convergenceTolerance,MAX_ITERS,mindx,forceBracketing)
    if ~exist('convergenceTolerance','var')
        convergenceTolerance = 1*10^-7;
    end
    if ~exist('forceBracketing','var')
        forceBracketing = false;
    end
    if ~exist('MAX_ITERS','var')
        MAX_ITERS = 1000;
    end
    if ~exist('mindx','var')
        mindx = 0.5*10^-16;
    end
    
    Upper = maxX;
    Lower = minX;

    yMin = func(minX);
    yMax = func(maxX);
    bracketing = true;
    if (yMin > 0 && yMax > 0) || (yMax < 0 && yMin < 0)
        if(forceBracketing)
            error('No change of sign within brackets and bracketing was forced!');
           %return; 
        end
        bracketing = false; %No sign change, cannot use bracketing
        x = fzero(func,x0); %Call to fzero instead
        return;
        %error('This solver is designed for when change of sign from xMin to xMax, problem supplied does not seem to satisfy this');
    end

    x = x0;
    y = func(x0); %Starting f(x)
    if bracketing
        dFdx = (yMax - yMin)/(maxX - minX); %d(func)/dx approximately - gradient of function solving
    else
        y1 = func(x0+convergenceTolerance);
        dFdx = (y1-y) / convergenceTolerance;
        %APPROX gradient at x0
    end
    dx = -((y) / dFdx); %Initial dx
    
    iter = 0; %Current iteration num
    while iter < MAX_ITERS
        if abs(dx) < mindx
            if dx < mindx
                dx = mindx;
            end
            if -dx < mindx
                dx = -mindx;
            end
        end
        iter = iter+1;
        yOld = y; %Keep a track of y from last iteration for approximating gradient of error function being solved
        x = x+dx; %Advance our guess for P

        %Ensure X is within bounds where change of sign exists
        if (x) < Lower %If newton raphson jumping around, try acting like interval bisection instead
            if bracketing
                dx = x-dx-Lower-(0.5*(Upper-Lower)); %Update dx to be the actual change in X from last iter
                x = Lower+(0.5*(Upper-Lower));
            else
                dx = x-dx-Lower;
                x = Lower;
            end
        end
        if (x) > Upper %If newton raphson jumping around, try acting like interval bisection instead
            if bracketing
                dx = Upper-(0.5*(Upper-Lower))-x-dx; %Update dx to be the actual change in X from last iter
                x = Upper-(0.5*(Upper-Lower));
            else
                dx = Upper-x-dx;
                x = Upper;
            end
        end

        y = func(x); %Calculate value of function at point
        %fprintf(['x: ',num2str(x),', y: ',num2str(y),'\n']);

        if dx == 0 %If dx=0,make non-zero so dFdX approximation is not Inf
            dx = 1; %y-yOld will also be 0
        end
        dFdx = (y-yOld)/dx; %Approximate rate of change of function
        if dFdx == 0 %If no approximation of gradient exists yet (Step too small?)
            %Use previous dx
        else
            dx = -((y) / dFdx); %Newton raphson x1=x0-f(x)/f'(x)
        end

        if abs(y) < convergenceTolerance %If close enough to solution
            return; %X found
        end

        if bracketing
            if x > Lower && ((y<0&&yMin<0) || (y>0&&yMin>0)) %If X is suitable as the new lower bound
                Lower = x;
                %fprintf(['Set ',num2str(x),' as lower limit\n']);
            elseif x < Upper && ((y<0&&yMax<0) || (y>0&&yMax>0)) %If X is suitable as the new upper bound
                %fprintf(['Set ',num2str(x),' as upper limit\n']);
                Upper = x;
            end
        end
        
%         disp("Lower: "+(Lower-184.6839)+" Upper: "+(Upper-184.6839));
    end

    %If reach here then didn't converge to a solution
    error('Unable to find required solution within max number of iterations');
end