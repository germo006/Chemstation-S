function [d1] = deriv1(X, z)
% This function will calculate first umerical derivatives to an
% error order of O(h^2)

    %INPUT%
    % X: matrix or vector of values where the derivatives are calculated
    % based on columns. This is not a 2D gradient, but several independent
    % variables.
    % z: either single value or vector of coordinates. Prefer even spacing.
    %
    %OUTPUT%
    % d1: matrix or vector, first derivatives using central difference or
    % three-point forward/backward differences for endpoints.
    
    
if length(z)==1
    h = ones(size(X,1)-1,1).*z; 
else
    h = diff(z);
end

d1 = zeros(size(X));

fDif = @(x0, x1, x2, h) (-3*x0 + 4*x1 - x2)./(2*h);
cDif = @(x1, xm1, h) (x1 - xm1)./(2.*h);
bDif = @(x0, xm1, xm2, h) (3*x0 - 4*xm1 + xm2)./(2*h);

d1(1,:) = fDif(X(1,:), X(2,:), X(3,:), mean(h(1:2)));
d1(2:end-1,:) = cDif(X(3:end,:), X(1:end-2,:), (h(1:end-1)+h(2:end))./2);
d1(end,:) = bDif(X(end,:), X(end-1,:), X(end-2,:), mean(h(end-1:end)));
