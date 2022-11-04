function [r2] = compute_rSquared(x,y)
%compute_rSquared computes the corfficient of determination for data x (observations) & y (predictions).
%For a least squares fit, y would be the predicted values
% Arrays are [samples x predictions]
% 
% Michael Rosenberg - U of Washington - 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preliminaries
assert(size(x,2) >= size(y,2)); % ensure Nx1 or 1xN vector
% If many predictions in x compared to one set of observations, repeat the
% observations
if size(x,2) > 1 && size(y,2) == 1 
    y = repmat(y,1,size(x,2));
end
    
N=length(x);
%% Compute r^2
xbar=sum(x,1)/N;
e=x-y;
if any(isnan(e))
    disp('Warning: Predictions or data contain NaN values that will be omitted from the r-squared calculation.')
end
SStot=nansum(( x - xbar ).^2,1)+eps; % Add a small amount of noise to avoid dividing by zero
% SSreg=nansum((y - xbar).^2, 1 );
SSres=nansum(e.^2, 1);
r2=1-SSres./SStot;
end

