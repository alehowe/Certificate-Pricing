function X=UpfrontLewisNTS(dates, discounts, ContractParams, NTSParams)
% UpfrontLewisNTS: computes the upfront of the contract relying on Lewis'
% closed formula

% Inputs:
% dates:            vector of the bootstrapped dates
% discounts:        vector of discounts at the corresponding dates
% ContractParams:   Struct containing the contract parameters 
% NTSParams:        Struct containing the NTS parameters

% Outputs:
% X:                Upfront of the contract (in bps)


% Define the parameters
Act360=2;
Act365=3;
Eu30_360=6;

% Extract the useful quantities
expiries=ContractParams.expiries;
EurPaymentDates = ContractParams.EurPaymentDates;
today=ContractParams.today;
strike=ContractParams.strike;
S0=ContractParams.S0;
q=ContractParams.q;
coupons=ContractParams.coupons;

sigma=NTSParams.sigma;
eta=NTSParams.eta;
k=NTSParams.k;
alpha=NTSParams.alpha;

% Compute the yearfracs between Euribor payment dates
aux_dates=[ContractParams.today, EurPaymentDates];
yearfracs_euribor=yearfrac(aux_dates(1:end-1), aux_dates(2:end), Act360);

% Compute the index up to the Euribor payment dates of the first year
% (reminding that expiries(1) contains the check date)
index_euribor=find(EurPaymentDates==expiries(2));

% Compute the discounts used at Euribor payment dates
used_discounts_euribor = find_discount(dates, discounts, EurPaymentDates);

% Compute the discounts at expiry
used_discounts_expiries=find_discount(dates, discounts, expiries);

% Compute the first zero rate
r=-log(used_discounts_expiries(1))./yearfrac(today, expiries(1),Act365);

% Compute the time to maturity up to the first check
ttm_check=yearfrac(today, expiries(1), Act365);

% Compute the N(d2) term: start with klewis
klewis=log(S0/strike)+(r-q)*ttm_check;

% Fuction of the Laplace transform
Lap_w= @(w) exp(ttm_check./k.*(1-alpha)./alpha.*(1-(1+(w.*k*sigma^2)./(1-alpha)).^alpha));

% Characteristic function
CharFct_psi=@(u) exp(-1i*u.*log(Lap_w(eta))).*Lap_w((u.^2+1i.*(1+2*eta).*u)/2);

% Find N(d2): the probability ot to trigger the value
tot_integrand = @ (u) real(exp(1i.*u.*klewis).* CharFct_psi(u)./(1i.*u));
N_d2 = 0.5+1/pi*quadgk(tot_integrand, 0, inf);
prob_first_coupon=1-N_d2;

% Compute the yearfracs used in the coupons
yearfracs_coupons=yearfrac([today, expiries(2)], [expiries(2), expiries(3)], Eu30_360);

% Compute the NPV of the coupon (reminding that the coupon are paid at the given dates)
NPVCoupon=prob_first_coupon*used_discounts_expiries(2)*yearfracs_coupons(1)*coupons(1)+(1-prob_first_coupon)*used_discounts_expiries(3)*coupons(2)*yearfracs_coupons(2);

% Compute the floating BPV 
BPVf_total=sum(yearfracs_euribor(1:end).*used_discounts_euribor(1:end));
BPVf_partial=sum(yearfracs_euribor(1:index_euribor).*used_discounts_euribor(1:index_euribor));

% Compute the NPV of Euribor 3m + Spol considering the trigger
NPVLib=(ContractParams.s_spol*BPVf_partial+1-used_discounts_euribor(index_euribor)).*prob_first_coupon+ ...
    (ContractParams.s_spol*BPVf_total+1-used_discounts_euribor(end)).*(1-prob_first_coupon);

% Compute the total NPV of the contract
X=(NPVLib-NPVCoupon)*1e4;

end