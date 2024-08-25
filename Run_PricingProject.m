% runAssignment7
% Howe Alessandro John

clear all;
close all;
clc;

%% Settings
formatData='dd/mm/yyyy'; 

%% Read market data
% This function works on Windows OS.

[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);


%% Bootstrap
% dates includes SettlementDate as first date

[dates, discounts]=bootstrap(datesSet, ratesSet); 

%% Calibrate the parameters

% Compute NIG parameters
alpha=0.5;
[sigma, eta, k]=computeNTSparameters(datesSet, dates, discounts, alpha);

% Store the parameters
NIGParams=struct();
NIGParams.sigma=sigma;
NIGParams.eta=eta;
NIGParams.k=k;
NIGParams.alpha=alpha;

% Generalization: compute the parameters in the generic normal tempered stable case (e.g. if
% alpha=1/3)
alpha=1/3;
[sigma, eta, k]=computeNTSparameters(datesSet, dates, discounts, alpha);

% Store the parameters
NTSParams=struct();
NTSParams.sigma=sigma;
NTSParams.eta=eta;
NTSParams.k=k;
NTSParams.alpha=alpha;

%% Storing contract's data

% Extract the relevant quantities
Act365=3;
dataset=load('cSelect20080215_B.mat');

% Build struct containing the contract parameters
ContractParams= struct();
ContractParams.today=dates(1);
ContractParams.maturities_changes=[1,2];
ContractParams.s_spol=0.013;
ContractParams.coupons=[0.06; 0.02];
ContractParams.X=[];
ContractParams.strike=3200;
ContractParams.notional=10^7;
ContractParams.S0=dataset.cSelect.reference;
ContractParams.q=dataset.cSelect.dividend;

% Build struct containing the MC parameters
MCParams=struct();
MCParams.M=1000000;
MCParams.seed=9;
MCParams.alpha=0.05;

% Build struct containing the volatility surface parameters
VolatilitiesParams=struct();
VolatilitiesParams.strikes=dataset.cSelect.strikes;
VolatilitiesParams.volatilities=dataset.cSelect.surface;

%% Pricing upfront

% Choose the flag:
% flag=1        for closed formula NIG model (using Lewis formula)
% flag=2        for MC simulation NIG model
% flag=3        for closed formula Black model (with digital risk)
% flag=4        for closed formula NTS model alpha=1/3 (using Lewis formula)
% flag=5        for MC simulation NTS model alpha=1/3 


% RMK. flag=4 and 5 can begeneralized to a generic alpha - you can pass the desired
% alpha in the function. However remember to recalibrate the parameters
% above and pass the corresponding struct

% For MC simulation (flag=2 and flag=5) also confidence interval (level alpha) and
% elapsed time is provided

flag=5;
% Compute X by setting the total NPV to 0
[X, CI_X, elapsed_time] = ComputeUpfrontContract(dates, discounts, ContractParams, NIGParams, ...
    NTSParams, MCParams, VolatilitiesParams, flag);
ContractParams.X=X;

% Print the contract's upfront
fprintf('---------------------Contract upfront------------------------\n');
fprintf('The contract''s upfront is: %.6f bps\n', X);
fprintf('Party B pays %.6f at start date\n', X/(10^4) * ContractParams.notional);


% Print confidence interval if Monte Carlo simulation
if ~isempty(CI_X)   
    fprintf('The confidence interval (in bps) for the upfront at level %.2f for %.0f Monte Carlo simulations is [%.6f, %.6f]\n', ...
        MCParams.alpha, MCParams.M, CI_X(1), CI_X(2));
    fprintf('The confidence interval for the upfront party B has to pay at level %.2f for %.0f Monte Carlo simulations is [%.6f, %.6f]\n', ...
        MCParams.alpha, MCParams.M, CI_X(1)/(10^4) * ContractParams.notional, CI_X(2)/(10^4) * ContractParams.notional);
end

% Print elapsed time if Monte Carlo simulation
if elapsed_time > 0  
    fprintf('Elapsed time for the MC simulation: %.6f seconds\n', elapsed_time);
end

%% Monte Carlo analysis

% Define a range of values for M
M_values = 10.^[3,4,5,6, 7, 8];
num_M_values = length(M_values);

% Initialize storage for results
X_values = zeros(num_M_values, 1);
CI_X_values = zeros(num_M_values, 2);
elapsed_times = zeros(num_M_values, 1);

% Fixed alpha value
alpha = MCParams.alpha;
flag=5;

% Loop over different values of M
for i = 1:num_M_values
    % Update MCParams with current M
    MCParams.M = M_values(i);
    
    % Compute upfront payment and confidence interval
    [X, CI_X, elapsed_time] = ComputeUpfrontContract(dates, discounts, ContractParams, NIGParams, ...
        NTSParams, MCParams, VolatilitiesParams, flag);
    
    % Store results
    X_values(i) = X;
    CI_X_values(i, :) = CI_X;
    elapsed_times(i) = elapsed_time;
end

% Plot the results
figure;
errorbar(M_values, X_values, X_values - CI_X_values(:, 1), CI_X_values(:, 2) - X_values, 'o-');
set(gca, 'XScale', 'log');
xlabel('Number of Monte Carlo Simulations (M)');
ylabel('Upfront Payment (bps)');
title('Upfront Payment and Confidence Intervals vs. Number of MC Simulations');
grid on;

figure;
plot(M_values, elapsed_times, 'o-');
set(gca, 'XScale', 'log');
xlabel('Number of Monte Carlo Simulations (M)');
ylabel('Elapsed Time (seconds)');
title('Elapsed Time vs. Number of MC Simulations');
grid on;

% Display results in a table
results_table = table(M_values', X_values, CI_X_values(:, 1), CI_X_values(:, 2), elapsed_times, ...
    'VariableNames', {'M', 'UpfrontPayment', 'CI_Lower', 'CI_Upper', 'ElapsedTime'});
disp(results_table);



