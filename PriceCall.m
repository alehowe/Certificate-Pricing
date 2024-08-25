function call_price= PriceCall(discount, alpha, sigma, nu,k, x, ttm, Fwd0,extraparam)
% PriceCall computes the price of a call option using the Lewis (2001) formula.

% If extraparam is not provided function computes the price of a call option
% at the given log-moneyness x using
% the Lewis (2001) formula using classical built-in quadrature MATLAB
% methods.

% If extraparam is provided function computes the price of a call 
% option with FFT integration using 
% the provided initialization parameters in extraparam.

% Inputs:
% discount:             Discount factor at expiry.
% alpha:                Alpha parameter of the normal tempered stable model.
% sigma:                Volatility parameter.
% nu:                   Nu parameter of the normal tempered stable model.
% k:                    Strike price.
% x:                    Log-moneyness at wich price is computed.
% ttm:                  Time to maturity.
% Fwd0:                 Initial forward price.
% extraparam:           Structure containing initialization parameters for FFT
%                       integration (optional).                  
%                       - For initialization with dz: extraparam contains
%                       fields M and dz.
%                       - For initialization with psi_start: extraparam
%                       contains fields M and psi_start.

% Outputs:
% call_price: Price of the call option.

% Function of the Laplace exponent in the normal tempered stable case
lnLap_w= @(w) ttm/k*(1-alpha)/alpha*(1-(1+(w.*k*sigma^2)/(1-alpha)).^alpha);

% Fuction of the Laplace transform
Lap_w= @(w) exp(ttm/k*(1-alpha)/alpha*(1-(1+(w.*k*sigma^2)/(1-alpha)).^alpha));

% Characteristic function
CharFct_psi=@(psi) exp(-1i*psi*lnLap_w(nu)).*Lap_w((psi+1i*(1+2*nu)).*psi/2);

% Integrand in Lewis' formula 
integrand= @(psi) 1/(2*pi)*CharFct_psi(-psi-1i/2)*1./(psi.^2+1/4);

% Compute the call price at a given log-moneyness 
if nargin()==9

    call_price=discount*Fwd0*(1-exp(-x/2).*solveintegral(x,integrand,extraparam));
    
else if nargin()==8

    call_price=discount*Fwd0*(1-exp(-x/2).*solveintegral(x,integrand));

else 
    error("Non valid integration method")
end

end