function X=UpfrontBlack(dates, discounts, ContractParams, VolatilitiesParams)
% UpfrontBlack: computes the upfront of the contract relying on adjusted
% Black's closed formula (i.e. considering digital risk)

% Inputs:
% dates:            vector of the bootstrapped dates
% discounts:        vector of discounts at the corresponding dates
% ContractParams:   Struct containing the contract parameters 
% VolatilitiesParams: Struct containing the volatility surface
%                     parameters

% Outputs:
% X:                Upfront of the contract


% Extract some quantities
expiries=ContractParams.expiries;
EurPaymentDates = ContractParams.EurPaymentDates;
today=ContractParams.today;
strike=ContractParams.strike;
S0=ContractParams.S0;
q=ContractParams.q;
coupons=ContractParams.coupons;

strikes=VolatilitiesParams.strikes;
volatilities=VolatilitiesParams.volatilities;

% Define the parameters
Act360=2;
Act365=3;
Eu30_360=6;

% Compute the yearfracs between Euribor payment dates
aux_dates=[ContractParams.today, EurPaymentDates];
yearfracs_euribor=yearfrac(aux_dates(1:end-1), aux_dates(2:end), Act360);

% Compute the index up to the Euribor payment dates of the first year
% (reminding that expiries(1) contains the check date)
index_euribor=find(EurPaymentDates==expiries(2));

% Compute the discounts used at Euribor payment dates
used_discounts_euribor = find_discount(dates, discounts, EurPaymentDates);

% Evaluate the discount factor 
discount_check=find_discount(dates,discounts,expiries(1)); 
discount_expiries=find_discount(dates,discounts,expiries(2:3));

% Compute F0 between today and the check step
ttm_codition_check = yearfrac(today,expiries(1),Act365);  % Act365!
r = -log(discount_check)./ttm_codition_check ;
F0_check = S0.*exp((r-q).*ttm_codition_check);

% Price the digital option considering the digital risk
smile_digital=DigitalPriceImpliedVol(F0_check, strike, strikes, volatilities, ttm_codition_check, discount_check, 1);

% Compute the probability of triggering the first coupon
prob_first_coupon=1-1/discount_check*smile_digital;

% Compute the yearfracs used in the coupons
yearfracs_coupons=yearfrac([today, expiries(2)], [expiries(2), expiries(3)] , Eu30_360);

% Compute the NPV (reminding that the coupon are paid at the given dates)
NPVCoupon=prob_first_coupon*discount_expiries(1)*yearfracs_coupons(1)*coupons(1)+(1-prob_first_coupon)*discount_expiries(2)*coupons(2)*yearfracs_coupons(2);

% Compute the floating BPV 
BPVf_total=sum(yearfracs_euribor(1:end).*used_discounts_euribor(1:end));
BPVf_partial=sum(yearfracs_euribor(1:index_euribor).*used_discounts_euribor(1:index_euribor));

% Compute the NPV of Euribor 3m + Spol considering the trigger
NPVLib=(ContractParams.s_spol*BPVf_partial+1-used_discounts_euribor(index_euribor)).*prob_first_coupon+ ...
    (ContractParams.s_spol*BPVf_total+1-used_discounts_euribor(end)).*(1-prob_first_coupon);

% Compute the total NPV of the contract
X=(NPVLib-NPVCoupon)*1e4;

end