%Non linear numerical equation solver for solving N non-linear equations
%equal to 0. Equations must take scalars as inputs eg. f(x,y,z) not f(x_vector)
%By Eddie Brown
classdef NonLinearEquationSolver < handle
    properties(Constant)
        SMALL_PERTUB = 1*10^-7; %Used for calculating approx partial derivatives
        DEFAULT_CONVERGENCE_TOLERANCE = 1*10^-6;
    end
    properties (SetAccess=private) %Not allowed to be externally modified
        equations = {};
    end
    properties %Allowed to be externally modified
        convergenceTolerance = NonLinearEquationSolver.DEFAULT_CONVERGENCE_TOLERANCE;
        MAX_ITERATIONS = 1000;
        MAX_STEP_SIZE_IN_ANY_DIMENSION = Inf;
    end
    properties (Dependent) %Dependent properties
        numEquations;
    end
    
    methods
        %Create a NonLinearEquationSolver object
        %Args: (equations,convergenceTolerance)
        %Equations is a cell array of function handles that comprise the
        %set of non-linear equations to solve
        %Convergence tolerance is how far from eqns all being 0 is allowed
        %for when considered solved.
        function obj = NonLinearEquationSolver(equations,convergenceTolerance) %Constructor for object
            obj.equations = equations;
            if ~exist('convergenceTolerance','var') || isempty(convergenceTolerance)
                obj.convergenceTolerance = NonLinearEquationSolver.DEFAULT_CONVERGENCE_TOLERANCE;
            else
                obj.convergenceTolerance = convergenceTolerance;
            end
        end
        
        %Gets function handle for equation number n contained by this
        %non-linear equation solver.
        function func=getEquation(obj,n)
            if n > obj.numEquations || n < 1
               error(['Provided equation number ',num2str(n),' invalid for number of equations this solver contains (',num2str(obj.numEquations),')']); 
            end
            
            func = obj.equations{n};
        end
        
        %Define getter function for number of equations
        function num=get.numEquations(obj)
            num = length(obj.equations);
        end
        
        %Attempt to solve the system of non-linear equations by starting
        %with the solution given for the given 'constants' and adjusting the
        %'constants' of the problem until they end up at the value we need
        %the solution for. 'constants' here mean any parameter not being
        %solved for. For example if you are solving for the Temp and
        %Velocity of a flow then your input 'constants' may be Po,To and P
        %and your known solution would be T=To,v=0 for when P=Po. This
        %function would then change Po in increments until it equals P,
        %using the solution at each step to predict a good starting point
        %guess for the next step. Use this function to solve problems where
        %convergence is only good for high quality guesses and you only know a
        %high quality guess somewhere else in the function than where your
        %problem is. stepInConstants is a function handle points to a
        %MATLAB function that takes the current value of the 'constants'
        %and the constants of the real problem eg. function(CCurrent,CReal)
        %and returns how much to move each one by so that it gets closer to
        %the solution.
        function X = solveFromKnownSolution(obj,knownSolution,constantsOfKnownSolution,constantsOfRealProblem,stepInConstants,lowerBound,upperBound)
            if ~exist('lowerBound','var') %if no lowerbound specified, set all to -Inf
               lowerBound = zeros(1,length(initialGuess)); 
               for p=1:length(lowerBound)
                  lowerBound(p) = -Inf; 
               end
            end
            if ~exist('upperBound','var') %if no lowerbound specified, set all to -Inf
               upperBound = zeros(1,length(initialGuess)); 
               for p=1:length(upperBound)
                  upperBound(p) = Inf; 
               end
            end
            
            constantsOfKnownSolution = obj.makeRowVector(constantsOfKnownSolution); %Ensure consistency
            constantsOfRealProblem = obj.makeRowVector(constantsOfRealProblem); %Ensure consistency
            
            %How much the solution changed by in the last 'step' (last
            %change of problem to be closer to the real one), divided by
            %the magnitude of the step
            dX = zeros(1,length(knownSolution));
            %Current value of 'constants'
            C = obj.makeRowVector(constantsOfKnownSolution);
            %Current solution
            X = obj.makeRowVector(knownSolution);
            currentStepSize = 1;
            while true
                isLast = C == constantsOfRealProblem; %True if the last iteration
                
                XOld = X;
                XGuess = XOld + dX*currentStepSize;
                attempt = 0;
                while(true) %Do this until it succeeds
                    try
                        %Try and solve
                        X = obj.solve(XGuess,C,lowerBound,upperBound);
                        break; %And stop looping
                    catch exception
                        %If solver failed, try adjusting initial guess
                        %before giving up
                        if strcmp(exception.identifier,'NonLinearEquationSolver:notConverge') || strcmp(exception.identifier,'NonLinearEquationSolver:failNaN') %Did not converge
                            attempt = attempt + 1;
                            
                            %Let's try moving a bit more and seeing if then
                            %converges (Some jumps in our data can cause it to get stuck in places)
                            XGuess = XGuess + dX*currentStepSize;

                            if attempt > 5 %If still not sorted, then just give up
                                disp("Max guess attempts (5) exceeded and still fails!");
                                disp("C:");
                                disp(C);
                                disp("X:");
                                disp(X);
                                rethrow(exception); %Rethrow error, stop execution
                            end
                        else
                            disp("Unhandled exception in solver! ("+exception.identifier+")");
                            rethrow(exception); %Rethrow error, stop execution
                        end
                    end
                end
                
                dX = (X - XOld)/currentStepSize;
                if all(dX==0) %If at initial solution, eg. X has not changed
                    dX(:) = 1*10^-5; %Give a little pertubation to our guess, helps
                end
                
                if isLast %If last iteration
                   return; %Finish
                end
                
                %Now decrease pressure by decrement amount
                currentStep = stepInConstants(C,constantsOfRealProblem);
                currentStepSize = norm(currentStep);
                C = C + currentStep;
                if dot(C-constantsOfRealProblem,C-constantsOfRealProblem) < currentStepSize^2 %If now have reached problem within distance of step
                   C = constantsOfRealProblem; %Calculate value at real problem next iteration
                end
            end
        end
        
        %Attempt to solve the system of non-linear equations with given
        %initialGuess. Note that with how this solver works, whether or not
        %the solver converges is largely determined by the quality of the
        %initial guess. See image in same folder as this for explanation of
        %how the principle behind this works. lowerBound is min values for
        %guess, upperBound is max values for guess. constants is a row
        %vector of constants to be passed to the functions evaluating each
        %equation, it is merged after the initialGuess vector. This is
        %useful is you need to pass parameters to the equations that are
        %not part of the problem (eg. they are constant). eg. initialGuess
        %[1,2] and constants [3,4,5] would pass to your function the inputs [1,2,3,4,5] but
        %would only change [1,2] in the solving process.
        function solution = solve(obj,initialGuess,constants,lowerBound,upperBound)
            warning('off','MATLAB:singularMatrix');
            warning('off','MATLAB:nearlySingularMatrix');
            if ~exist('constants','var')
                constants = [];
            else
               constants = obj.makeRowVector(constants); 
            end
            if ~exist('lowerBound','var') %if no lowerbound specified, set all to -Inf
               lowerBound = zeros(1,length(initialGuess)); 
               for p=1:length(lowerBound)
                  lowerBound(p) = -Inf; 
               end
            end
            if ~exist('upperBound','var') %if no lowerbound specified, set all to -Inf
               upperBound = zeros(1,length(initialGuess)); 
               for p=1:length(upperBound)
                  upperBound(p) = Inf; 
               end
            end
            
            if length(initialGuess) < obj.numEquations
               error('Initial guess supplied to equation solver did not have enough unknowns to be solvable'); 
            end
            x0 = obj.makeRowVector(initialGuess);
            x = x0;
            iterNum = 0;
            xOld = x;
            while (true)
                %Enforce that 'guess' is within allowed limits of domain
               for k=1:length(x)
                  if x(k) > upperBound(k) %If too large
                      x(k) = lowerBound(k); %Likely wasn't converging, so jump to other side of domain in hope it'll converge
                  elseif x(k) < lowerBound(k) %If too small
                      x(k) = upperBound(k); %Likely wasn't converging, so jump to other side of domain in hope it'll converge
                  end
               end
               
               iterNum = iterNum+1;
               if iterNum > obj.MAX_ITERATIONS
                   throw(MException('NonLinearEquationSolver:notConverge','Max number of iterations reached and solution did not converge!'));
               end
               
               if any(isnan(x))
                   throw(MException('NonLinearEquationSolver:failNaN',['Solver failed to solve, ended up with NaN values in x in iteration ',num2str(iterNum),'!']));
               end
               
               f = obj.evalAllEquations([x,constants]); %Current val of every equation with current guess for x, as a row vector
               if all(abs(f) < obj.convergenceTolerance) %If ALL f(i) values are within convergence tolerance of zero
                   solution = x; %Consider the system solved! WOOHOO!
                   warning('on','MATLAB:singularMatrix');
                   warning('on','MATLAB:nearlySingularMatrix');
                   return;
               end
               
               %Solve linear system J*h = -f where J is Jacobian matrix, h is col
               %vector of how much to move each x coordinate by and f is a
               %col vector of the current function evaluation
               J = obj.computeJacobian(x,constants);
               h = linsolve(J,-transpose(f));
               
               for r=1:length(h)
                  if(abs(h(r)) > obj.MAX_STEP_SIZE_IN_ANY_DIMENSION)
                      if h(r) < 0
                          h(r) = -obj.MAX_STEP_SIZE_IN_ANY_DIMENSION;
                      else
                          h(r) = obj.MAX_STEP_SIZE_IN_ANY_DIMENSION;
                      end
                  end
               end
               
               xOldOld = xOld;
               xOld = x;
               x = x+transpose(h); %Improve guess for x
               
               if any(isnan(f)) %If evaluation of functions results in NaNs
                   for l=1:length(x) %Revert to non-NaN x, with a bit of random deviation
                       xOld = xOldOld;
                      x(l) = xOldOld(l)+((rand()-0.5)/5);
                    end
               end
               for l=1:length(x)
                  if isnan(x(l)) || isinf(x(l)) %Nan or inf val in this x, gradient at this point probably a problem
                      x(l) = xOld(l)+((rand()-0.5)/5); %Hopefully small pertubation to guess will fix it
                  end
               end
            end
        end
        
        %Evaluates all equations for the given input vector and returns as
        %a row vector
        function yVals=evalAllEquations(obj,inputs)
            yVals = zeros(1,obj.numEquations); %Row vector of zeros
            for i=1:obj.numEquations
               yVals(i) = evalEquation(obj,i,inputs); %Eval yi and put into y vector 
            end
        end
        
        %Evaluates given equation number (n) with the provided inputs and
        %returns the result. inputs should be a vector with the inputs of
        %the function, eg. [x,y,z]
        function y=evalEquation(obj,n,inputs)
            func = obj.getEquation(n); %Get eqn
            inp = obj.toRowCellArr(inputs); %Convert input to required format for matlab to eval
            y = func(inp{:}); %Eval function with these inputs
        end
    end
    methods(Access=private)
        
        %Compute the jacobian matrix for the set of equations at the
        %position specified by x0 vector where x0 vector = [x0,x1,x2,...xN]
        %Jacobian matrix: https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
        function jacob=computeJacobian(obj,x0,constants)
            x0 = obj.makeRowVector(x0); %Ensure orientation is consistent
            jacob = zeros(obj.numEquations,length(x0)); %Jacobian matrix to populate
            for eqnNum=1:obj.numEquations %For each eqn (fi)
               for dimen=1:length(x0) %For each variable (xj)
                   jacob(eqnNum,dimen) = obj.approxGradInDimenOfFunc(eqnNum,[x0,constants],dimen);
               end
            end
        end
        
        %Approximates the partial derivative dfi/dxj for i=eqnNum and
        %j=dimenNum at the position given by the vector 'inputs' of the form [x0,x1,x2,...,xN]
        function grad=approxGradInDimenOfFunc(obj,eqnNum,inputs,dimenNum)
            inputs = obj.makeRowVector(inputs); %Ensure orientation is consistent
            dxi = NonLinearEquationSolver.SMALL_PERTUB;
            y1 = obj.evalEquation(eqnNum,inputs);
            inputs(dimenNum) = inputs(dimenNum) + dxi;
            y2 = obj.evalEquation(eqnNum,inputs);
            grad = (y2-y1) / dxi; %Approximate gradient of function
        end
        
        function vec = makeRowVector(obj,v)
            sizeInp = size(v);
            numVars = length(v);
            if sizeInp(1) >  numVars %If is a column vector
                v = transpose(v); %Make it a row vector, for consistency
            end
            vec = v;
        end
        
        %Function that converts from matlab vector to a row cell array, so can be
        %passed to a function handle for example
        function cellArr = toRowCellArr(obj,x0)
            %Convert x0 to cell array
            sizeInp = size(x0);
            numVars = length(x0);
            if sizeInp(1) >  numVars %If is a column vector
                x0 = transpose(x0); %Make it a row vector, for consistency
                sizeInp = size(x0);
            end
            
            xWrapped = {};
            for i=1:length(x0)
                xWrapped(i) = mat2cell(x0(i), 1, 1);
            end
            
            cellArr = xWrapped;
        end
    end
    
end