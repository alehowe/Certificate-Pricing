function [X, CI_X, elapsed_time] =ComputeUpfrontContract(dates, discounts, ContractParams, NIGParams, NTSParams, MCParams, VolatilitiesParams, flag)
% ComputeUpfrontContract: computes the upfront of the certificate

% Input:
% dates:          Vector of dates.
% discounts:      Vector of discount factors corresponding to the dates.
% ContractParams: Struct containing the contract parameters 
% NIGParams:      Struct containing the NIG parameters
% NTSParams:      Struct containing the NTS parameters
% MCParams:       Struct containing the MC paramters
% VolatilitiesParams: Struct containing the implied surface
%                     volatilities parameters
% flag:           Flag which determines the coupon pricing method:    
%                 1: Closed formula using NIG - Lewis formula
%                 2: NIG MC simulation 
%                 3: Closed formula using Black (with digital risk)
%                 4: Closed formula using NTS for given alpha - Lewis formula
%                 5: MC simulation using NTS for given alpha
 
% Outputs:
% X:                Upfront of the contract in bps
% CI_X:             (ONLY IF MONTE CARLO): Confidence interval at level 
%                   alpha (specified in MCParams) of the upfront in bps
%                   
% elapsed_time:     (ONLY IF MONTE CARLO): Elapsed time for MC simulation 
%                   in seconds


% Initialize CI_X and elapsed_time
CI_X= [];
elapsed_time = -1;

% Parameters
Act365=3;
Act360=2;
Eu_30_360=6;

% Initialize the used quantities
EurPaymentDates=zeros(1, ContractParams.maturities_changes(end)*4);
expiries=zeros(1, ContractParams.maturities_changes(end));

% Compute the value dates starting from today up to the first year
for  ii=1:ContractParams.maturities_changes(end)*4

    % Dates every 3 months from now to expiry
    EurPaymentDates(ii)=addtodate(ContractParams.today, 3*ii, 'month');
end

for  ii=1:ContractParams.maturities_changes(end)

    % Dates every 3 months from now to expiry
    expiries(ii)=addtodate(ContractParams.today, ii, 'year');
end


% Consider only the business dates
%holidays=holidays(today, EurPaymentDates(end)+1);
EurPaymentDates=busdate(EurPaymentDates-1, 'follow', holidays);


% Add the 2 days before the maturity to the expiries and check it isn't a
% hoiday
dates_2_days_before=addtodate(expiries(1), -2, 'day');
while ~isbusday(dates_2_days_before)
    dates_2_days_before=addtodate(dates_2_days_before, -1, "day");
end

expiries=sort([expiries, dates_2_days_before]);

% Store these variables in the struct
ContractParams.expiries=expiries;
ContractParams.EurPaymentDates=EurPaymentDates;

% Compute the Coupon NPV according to the chosen flag
switch(flag)
    case 1
        X = UpfrontLewisNTS(dates, discounts, ContractParams, NIGParams);
    case 2
        [X, CI_X, elapsed_time]=UpfrontMCNTS(dates, discounts, ContractParams, NIGParams, MCParams);
    case 3 
         X=UpfrontBlack(dates, discounts, ContractParams, VolatilitiesParams);
    case 4
         X=UpfrontLewisNTS(dates, discounts, ContractParams, NTSParams);
    case 5
        [X, CI_X, elapsed_time]=UpfrontMCNTS(dates, discounts, ContractParams, NTSParams, MCParams);
    otherwise
        error("Flag not accepted")
end


end