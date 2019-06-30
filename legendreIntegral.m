%Function to calculate the legendre integral of a function for a given
%number of legendre quadrature points with given limits a and b
function I = legendreIntegral(func,numPts,a,b)
    %Use change of variables x=0.5(a+b) + 0.5(b-a)t, dx=0.5(b-a)dt
    x = @(t) 0.5.*(a+b) + 0.5.*(b-a).*t;
    f2 = @(t) func(x(t)) .* 0.5 .* (b-a); %Integrand of gauss legendre
    [X,W] = legpts(numPts); %Find gauss legendre points and weights via chebfun script off the internet
    I = 0;
    for i=1:numPts
        I = I+W(i)*f2(X(i)); %Weight * integrand at point
    end
end