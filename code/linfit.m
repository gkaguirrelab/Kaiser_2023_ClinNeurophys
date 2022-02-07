function [Coeffs, Rho, R2, Residuals] = linfit(x,y)
% call:
% 
%     [Coeffs, Rho, R2, Residuals] = linfit(x,y)
%     
% Compute linear fitting between vectors 'x' and 'y', and return the coefficients of the regression 
% line, the correlation between 'x' and 'y', the coefficient of determination ("how much knowing x it
% is possible to determine y"), and the residuals. If the residuals exhibit "non-random" patterns, it 
% is an indication that a linear fit is not so good. To check them: bar(x,Residuals), xlabel('x'), ylabel('ŷ - y') 
% Note: the function ignores non-finite values(-Inf,Inf,NaN).

% INPUT
% 
%     x          : independent variable
%     y          : dependent variable
%     
% OUTPUT
% 
%     Coeffs     : Coefficients of the best-fitting linear regression line
%     Rho        : Correlation between 'x' and 'y'
%     R2         : Coefficient of Determination
%     Residuals  : difference between observations predicted by the best linear fit and the real ones
%     
%     
%
%
% Ruggero G. Bettinardi (RGB)
% Cellular & System Neurobiology, CRG
% -------------------------------------------------------------------------------------------

% Code History:
% --------------------------------------
%  ...  2015, RGB: Function is created
% 3 Feb 2017, RGB: Function is updated



% evaluate inputs:
% -----------------
if ~isvector(x) || ~isvector(y) 
    error('x1 and x2 should be vectors !');
end

nx = length(x);
ny = length(y);

if nx ~= ny
    error('x and y should have same length!');
end

% make the two vector have same size:
if size(x) ~= size(y),
    x = reshape(x,nx,1);
    y = reshape(y,nx,1);
end


if numel(x(isfinite(x))) ~= numel(x) || numel(y(isfinite(y))) ~= numel(y),
    display(' --->  NOTE: Non-finite values (e.g. Inf, -Inf, NaN) have been ignored')
    xfid  = find(isfinite(x)==1);  % indices of finite values in x
    yfid  = find(isfinite(y)==1);  % indices of finite values in x
    xyfid = intersect(xfid,yfid);  % indices of finite elements shared by x and y
    x     = x(xyfid);              % new vector x containing only shared non-finite values
    y     = y(xyfid);              % new vector x containing only shared non-finite values   
end
    

Rho       = corr(x,y);
Coeffs    = polyfit(x,y,1);
yhat      = polyval(Coeffs,x); 
Residuals = yhat - y;
SSR       = sum(Residuals.^2);
SST       = sum((y - mean(y)).^2); 
R2        = 1 - (SSR/SST);
