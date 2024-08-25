function [c, ceq] = NonlinearConstraintEta(params, alpha)
% NonlinearConstraintEta: defines the nonlinear constraint function to 
% impose on eta in the normal tempered stable case

% Inputs:
% params:           variable of the function (thus params=(sigma, eta, k))
% alpha:            alpha parameter in the normal tempered stable case

% Outputs:
% c:                inequality constraint
% ceq:              equality constraint


% Impose that eta>=-w_bar       (params=(sigma, eta, k))
c=-(1-alpha)/(params(3)*params(1)^2)-params(2);

% No equality to declare
ceq=[];

end
