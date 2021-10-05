function [cf,G]=L3P(x,y,varargin)
%L3P Three Parameters logistic regression
% The Three Parameters Logistic Regression or 3PL nonlinear regression model
% is commonly used for curve-fitting analysis in bioassays or immunoassays
% such as ELISAs or dose-response curves. 
% It is characterized by it’s classic “S” or sigmoidal shape that fits the
% top plateaus of the curve, the EC50, and the slope factor
% (Hill slope). This curve is symmetrical around its inflection point. 
%
% The 3PL equation is:
% F(x) = D/(1+(x/C)^(-B))
% where:
% B = Hill's slope. The Hill's slope refers to the steepness of the curve.
% It could either be positive or negative.
%
% C = Inflection point. The inflection point is defined as the point on the
% curve where the curvature changes direction or signs. C is the
% concentration of analyte where y=D/2.
%
% D = Maximum asymptote. In an bioassay where you have a standard curve,
% this can be thought of as the response value for infinite standard
% concentration. 
% 
%
% Syntax: [cf,G]=L3P(x,y,st,L,U)
% 
% Inputs: 
%           X and Y (mandatory) - data points.
%           X is a Nx1 vector and Y must have the same rows number
%           of X. If Y is a NxM matrix (N points and M replicates), L3P
%           will generate a  column vector computing means for each row.
%           The standard deviations of the rows will be used as weights of
%           regression.
%
%           st = starting points. This is a 1x3 vector of starting points
%           that have to be used to start the process of not linear
%           fitting. If this vector is not provided, L3P will set the
%           starting points on the basis of x and y data.
%
%           L = Lower bounds of parameters. This is a 1x3 vector of lower
%           bounds of the 3 parameters. If this vector is not provided, L3P
%           will set it on the basis of x and y data.
%
%           U = Upper bounds of parameters. This is a 1x3 vector of upper
%           bounds of the 3 parameters. If this vector is not provided, L3P
%           will set it on the basis of x and y data.
%
% Outputs:
%           cf = the FIT object
%           G = goodness-of-fit measures, for the given inputs, in the
%           structure G. G includes the fields: 
%           -- SSE         sum of squares due to error
%           -- R2          coefficient of determination or R^2
%           -- adjustedR2  degree of freedom adjusted R^2
%           -- stdError    fit standard error or root mean square error
% 
% Example:
%
% x=[0 4.5 10.6 19.7 40 84 210]'; y=[0.0089 0.0419 0.0873 0.2599 0.7074 1.528 2.7739]';
%
% Calling on MatLab the function: [cf,G]=L3P(x,y)
% 
%           Answer is:
% 
% cf = 
% 
%      General model:
%      cf(x) = D/(1+(x/C)^(-B))
%      Coefficients (with 95% confidence bounds):
%        B =       1.512  (1.384, 1.641)
%        C =       108.2  (93.29, 123.1)
%        D =       3.789  (3.471, 4.107)
% 
% G = 
% 
%            sse: 0.0012
%        rsquare: 0.9998
%            dfe: 4
%     adjrsquare: 0.9997
%           rmse: 0.0173
%
% hold on; plot(x,y,'ro'); plot(cf,'r'); hold off
% this will plot the curve.
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% See also L3Pinv, L4P, L4Pinv, L5P, L5Pinv
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2012) Three parameters logistic regression - There and back again
% https://it.mathworks.com/matlabcentral/fileexchange/38124

%--------------------Input errors handling section-------------------------
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'column','real','finite','nonnan','nonempty'}));
addRequired(p,'y',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty'}));
addOptional(p,'st',[],@(x) validateattributes(x,{'numeric'},{'row','real','finite','nonnan','ncols',3}));
addOptional(p,'L',[],@(x) validateattributes(x,{'numeric'},{'row','real','nonnan','ncols',3}));
addOptional(p,'U',[],@(x) validateattributes(x,{'numeric'},{'row','real','nonnan','ncols',3}));
parse(p,x,y,varargin{:});
st=p.Results.st; L=p.Results.L; U=p.Results.U;
clear p
assert(size(x,1)==size(y,1),'X and Y must have the same rows number')
assert(size(x,1)>=3,'Not enough Data points')

% if y is a matrix, compute means and standard deviations
if ~isvector(y) 
    we=std(y,0,2);
    y=mean(y,2);
else
    we=zeros(size(x));
end

%set the starting points:
% B is the Hill's slope so guess it with the slope of the line between first and last point.
% C is the inflection point (the concentration of analyte where you have
% half of the max response) so guess it finding the concentration whose
% response is nearest to the mid response.
% D is the upper asymptote so guess it with max(y)

slope=(y(end)-y(1))/(x(end)-x(1));
if isempty(st)
    [~,Idx]=min(abs((y-max(y)/2)));
    st=[sign(slope) x(Idx) max(y)];
end

%set the bounds. Of course all lower bounds are 0 and all upper bound are
%Inf. Anyway, if the slope is negative the lower bound of B is -Inf and the
%upper bound is 0. 

if isempty(L)
    L=zeros(1,3);
    if slope<0
        L(1)=-Inf;
    end
end

if isempty(U)
    U=Inf(1,3);
    if slope<0
        U(1)=0;
    end
end

%-----------------------------Fit the data---------------------------------
fo = fitoptions('method','NonlinearLeastSquares','Lower',L,'Upper',U);
set(fo,'Startpoint',st);
if all(we) % if y was a matrix use std as weights for fitting
    set(fo,'Weights',we);
end
ft = fittype('D/(1+(x/C)^(-B))',...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'B', 'C', 'D'});
 [cf,G] = fit(x,y,ft,fo);
