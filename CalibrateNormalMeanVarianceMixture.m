function [sigma_cal, eta_cal, k_cal]=CalibrateNormalMeanVarianceMixture(F0, strikes, surface, TTM, expirydiscount, alpha, weights, initial_guess, extra_params)
% CalibrateNormalMeanVarianceMixture: global calibration of the parameters 
% sigma, eta, k modelling the log forward dynamics as a normal variance mixture 
% given the implied volatility surface

% Inputs:
% F0:                   initial forward
% strikes:              vector containing different strikes
% surface:              vector containing the volatility surface corresponding
%                       to the vector K
% TTM:                  time to maturity of the option
% expirydiscount:       discount at expiry (considering TTM)
% alpha:                alpha parameter in the normal tempered stable case
% weights:              weights used for the global calibration
% initial_condition:    initial minimization point
% extra_params:         struct containing extra parameters in case of use of 
%                       the fft algorithm

% Outputs:
% sigma_cal:            estimated sigma
% eta_cal:              estimated eta
% k_cal:                estimated k


% Initialize the market prices
mkt_prices=zeros(1,length(strikes));

% Compute the zero rate
r=-log(expirydiscount)/TTM;

% Compute the market prices
mkt_prices=blkprice(F0, strikes, r, TTM, surface);

% Compute the log moneyness
x=log(F0./strikes);

% Compute the call price depending on the method
%params=[sigma, eta, k]; 
switch(nargin)
    case 8        
        normal_mean_mixture_price = @(params1,params2,params3) PriceCall(expirydiscount, alpha, params1, params2, params3, x, TTM, F0);
        
    case 9

        normal_mean_mixture_price = @(params1,params2,params3) PriceCall(expirydiscount, alpha, params1, params2, params3, x, TTM, F0, extra_params);

    otherwise
        error("Number of parameters don't match");
end

%sigma, eta,k,

% Define the objective function
to_minimize = @(params) sum(weights .* (mkt_prices - normal_mean_mixture_price(params(1), params(2), params(3))).^2);


% Define the linear lower bounds for the parameters
lb_sigma = 0;
lb_k = 0;

% Perform constrained optimization using fmincon (using
% NonlinearConstraintEta as the anonymous function)
% optimal_params = fmincon(to_minimize, initial_guess, [], [], [], [], [lb_sigma, -Inf, lb_k], [Inf, Inf, Inf], @NonlinearConstraintEta);
optimal_params = fmincon(@(params) to_minimize(params), initial_guess, [], [], [], [],...
    [lb_sigma, -Inf, lb_k], [Inf, Inf, Inf],@(params) NonlinearConstraintEta(params, alpha));
sigma_cal=optimal_params(1);
eta_cal=optimal_params(2);
k_cal=optimal_params(3);

end

