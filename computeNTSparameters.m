function [sigma_cal_quad, eta_cal_quad, k_cal_quad]=computeNTSparameters(datesSet, dates, discounts, alpha)
% computeNTSparameters: calibrates the parameters sigma, eta, k of a
% generic normal tempered stable model for a given alpha

% Inputs:
% datesSet:         dates struct
% dates:            vector containing the bootstrap dates
% discounts:        vector containing the discounts from the bootstrap
%                   corresponding to the given dates
% alpha:            alpha value of the normal tempered stable case

% Outputs:
% sigma:            sigma parameter
% eta:              eta parameter
% k:                k parameter
clc;
fprintf("Parameters calibration -------------------------------------------------\n")

% Set year fraction conventions
Act365=3;

% Load the smile
dataset=load('cSelect20080215_B.mat');

% Set today's date
today=datesSet.settlement;

% Compute expiry date (it is not said to consider business dates only)
expiry=addtodate(today, 1, 'year');

% Compute time to maturity
TTM=yearfrac(today,expiry,Act365);

% Evaluate the discount factor 
expirydiscount=find_discount(dates,discounts,expiry);

% Retrieve values
S0=dataset.cSelect.reference;
d=dataset.cSelect.dividend;
strikes=dataset.cSelect.strikes;
surface=dataset.cSelect.surface;

% Compute F0
r=-log(expirydiscount)/TTM;
F0=S0*exp(TTM*(r-d));

% Compute the log moneyness
x=log(F0./strikes);

% Define needed parameters (I use unitary weights)
weights=ones(1, length(strikes));  
initial_guess=[0.2, 3, 1];       %(sigma, eta, k) 

%% Calibrate using quadrature
% Calibrate sigma, eta, k according to quadrature method
[sigma_cal_quad, eta_cal_quad, k_cal_quad]=CalibrateNormalMeanVarianceMixture(F0,strikes, surface, TTM, expirydiscount, alpha, weights, initial_guess);

fprintf(['The calibrated parameters according to Quadrature integral computation are: \n ' ...
    'sigma: %f \n eta: %f \n, k: %f'], sigma_cal_quad, eta_cal_quad, k_cal_quad);

% Compute the implied volatility surface of the calibrated model
normal_mean_mixture_price_quad=PriceCall(expirydiscount, alpha, sigma_cal_quad, eta_cal_quad, k_cal_quad, x, TTM, F0);

% Plot the market and the model0's implied volatility surfaces 
PlotImpliedVolSurfaceCompared(F0, strikes, surface, TTM, expirydiscount, normal_mean_mixture_price_quad)
format long
% Compute the normalized errors
norm_errors=abs(blkimpv(F0, strikes, r, TTM, normal_mean_mixture_price_quad)-surface)./surface;

% Compute the strikes at which the error is maximum
[top_norm_errs, sorted_indices] = sort(norm_errors, 'descend');
top_5_errors = strikes(sorted_indices(1:5));

% Compute the strikes at which the error is minimum
[least_norm_errs, sorted_indices] = sort(norm_errors);
least_5_errors = strikes(sorted_indices(1:5));

% Compute the mean error
mean_error_quad=mean(norm_errors);
fprintf(['\nThe mean error in the implied volatility surface according to ' ...
    'Quadrature integral computation is: %f'], mean_error_quad)


