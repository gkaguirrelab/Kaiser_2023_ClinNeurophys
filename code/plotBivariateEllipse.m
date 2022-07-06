function plotBivariateEllipse(xVals,yVals,P)

% Observed
xMean = mean(xVals); % T mean
yMean= mean(yVals); % Td mean
CV = cov(xVals,yVals); % covariance of T and Td
[Evec, Eval]=eig(CV); % Eigen values and vectors of covariance matrix

%%% Plot observed multivariate contours %%

% Observed data
xCenter = xMean; % ellipses centered at sample averages
yCenter = yMean;
theta = 0 : 0.01 : 2*pi; % angles used for plotting ellipses

% compute angle for rotation of ellipse
% rotation angle will be angle between x axis and first eigenvector
x_vec= [1 ; 0]; % vector along x-axis
cosrotation =dot(x_vec,Evec(:,1))/(norm(x_vec)*norm(Evec(:,1))); 
rotation =pi/2-acos(cosrotation); % rotation angle
R  = [sin(rotation) cos(rotation); ...
      -cos(rotation)  sin(rotation)]; % create a rotation matrix

% create chi squared vector
chisq = chi2inv(P,2); % percentiles of chi^2 dist df=2

% size ellipses for each quantile
for i = 1:length(chisq)
    % calculate the radius of the ellipse
    xRadius(i)=(chisq(i)*Eval(1,1))^.5; % primary
    yRadius(i)=(chisq(i)*Eval(2,2))^.5; % secondary
    % lines for plotting ellipse
    x{i} = xRadius(i)* cos(theta);
    y{i} = yRadius(i) * sin(theta);
    % rotate ellipse
    rotated_Coords{i} = R*[x{i} ; y{i}];
    % center ellipse
    x_plot{i}=rotated_Coords{i}(1,:)'+xCenter;
    y_plot{i}=rotated_Coords{i}(2,:)'+yCenter;
end

% Plot contours
for j = 1:length(chisq)
    plot(x_plot{j},y_plot{j},'-k')
end

end