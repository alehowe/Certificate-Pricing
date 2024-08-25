function [X, CI_X, elapsed_time]=UpfrontMCNTS(dates, discounts, ContractParams, NTSParams, MCParams)
% UpfrontMCNTS: computes the upfront of the contract using a Monte Carlo
% simulation of the forward dynamics (whose exponent is modelled as a
% normal mean variance mixture, with the Laplace exponent of the its 
% characteristic funcrion is of a tempered stable positive random variable)

% Inputs:
% dates:            vector of the bootstrapped dates
% discounts:        vector of discounts at the correspondning dates
% ContractParams:   Struct containing the contract parameters 
% NTSParams:        Struct containing the NTS parameters
% MCParams:         Struct containing the MC parameters

% Outputs:
% X:                Upfront of the contract in bps
% CI_X:             Confidence interval of level alpha (specified in
%                   MCParams) of the upfront in bps
% elapsed_time:     Elapsed time for MC simulation in seconds


% Define the parameters
Act360=2;
Act365=3;
Eu30_360=6;

% Extract used quantities
expiries=ContractParams.expiries;
EurPaymentDates = ContractParams.EurPaymentDates;
today=ContractParams.today;
strike=ContractParams.strike;
S0=ContractParams.S0;
q=ContractParams.q;
coupons=ContractParams.coupons;
s_spol=ContractParams.s_spol;

sigma=NTSParams.sigma;
eta=NTSParams.eta;
k=NTSParams.k;
alpha=NTSParams.alpha;

M=MCParams.M;
seed=MCParams.seed;
alpha=MCParams.alpha;

% Compute the yearfracs between Euribor payment dates
aux_dates=[ContractParams.today, EurPaymentDates];
yearfracs_euribor=yearfrac(aux_dates(1:end-1), aux_dates(2:end), Act360);

% Compute the index up to the Euribor payment dates of the first year
% (reminding that expiries(1) contains the check date)
index_euribor=find(EurPaymentDates==expiries(2));

% Compute the discounts used at Euribor payment dates
used_discounts_euribor = find_discount(dates, discounts, EurPaymentDates);

% Compute the discounts at payment dates
used_discounts_coupons=find_discount(dates, discounts, expiries(2:3)); % at 1st and 2nd year
discount_check=find_discount(dates, discounts, expiries(1)); % at check date

% Start timing MC simulation
tic;

% Set the seed 
rng(seed);

% Simulate M standard normal random variables
g=randn(M,1);

% Simulate G
u=rand(M,1);
z=chi2rnd(1,M,1);
G=1-k/2*(sqrt(z.^2+4*z./k)-z);
G(find((1+G).*u>1))=1./G(find((1+G).*u>1));

% Compute F0 between today and the check step
ttm_codition_check = yearfrac(today,expiries(1),Act365);  % Act365!
r = -log(discount_check)./ttm_codition_check ;
F0_check = S0.*exp((r-q).*ttm_codition_check);

% Simulate the forward at the check step 
ln_lap_w= ttm_codition_check./k.*(1-alpha)./alpha.*(1-(1+eta.*k*sigma^2./(1-alpha)).^alpha);
f=sqrt(ttm_codition_check)*sigma.*sqrt(G).*g-(0.5+eta)*ttm_codition_check*sigma^2.*G-ln_lap_w;
F=F0_check.*exp(f);

% Stop timing after having simulated
elapsed_time = toc;

% Compute the yearfracs used in the coupons
yearfracs_coupons=yearfrac([today, expiries(2)], [expiries(2), expiries(3)], Eu30_360);

% Compute confidence interval of level alpha for the NPV of the coupons
C=ones(M,1);
C(F<strike,1)=coupons(1)*used_discounts_coupons(1)*yearfracs_coupons(1);
C(F>=strike,1)=coupons(2)*used_discounts_coupons(2)*yearfracs_coupons(2);

% Compute the floating BPV 
BPVf_total=sum(yearfracs_euribor(1:end).*used_discounts_euribor(1:end));
BPVf_partial=sum(yearfracs_euribor(1:index_euribor).*used_discounts_euribor(1:index_euribor));

% Compute the NPV for the floating leg in each scenario
EUR = ones(M,1);
EUR(F>=strike, 1) = s_spol*BPVf_total+1-used_discounts_euribor(end);
EUR(F<strike, 1) = s_spol*BPVf_partial+1-used_discounts_euribor(index_euribor);

% Compute the upfront in each scenario
UP = EUR-C;

% Compute the mean upfront
X=sum(UP)/M;

% Provide the confidence interval of level alpha for the upfront
sample_variance = sum((UP-X).^2)/(M-1);
CI_X = [-norminv(1-alpha/2)*sqrt(sample_variance/M)+ X; norminv(1-alpha/2)*sqrt(sample_variance/M)+ X];

% Return in bps
X = X*1e4;
CI_X = CI_X*1e4;

end