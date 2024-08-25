function digital_price=DigitalPriceImpliedVol(F0, K, strikes, volatilities, TTM, expirydiscount, notional)
% DigitalPriceImpliedVol: computes price of a digital option according to
% the implied volatility approach

% Inputs:
% F0:               initial forward
% K:                strike of the option
% strikes:          vector containing different strikes
% volatilities:     vector containing the volatility surface corresponding
%                   to the vector K
% TTM:              time to maturity of the option
% expirydiscount:   discount at expiry (considering TTM)
% notional:         notional 

% Outputs:
% digital_price:    price of the digital option


% Compute the slope impact passing by the two closest points to K in the dataset 
min_index=find(strikes>K, 1);
max_index=find(strikes>=K, 1)-1;

if isempty(min_index) || isempty(max_index)
    error('The strike is out of bounds');
end

slope_imp=(volatilities(max_index)-volatilities(min_index))/(strikes(max_index)-strikes(min_index));

% Compute the volatility correspondig to the given strike
digital_vol=interp1(strikes,volatilities,K,'spline');

% Compute the vega
d1 = 1./sqrt(TTM.*digital_vol^2).*(log(F0/K) + 1/2*TTM*digital_vol^2);
vega=expirydiscount*F0*exp(-d1.^2/2)/sqrt(2*pi)*sqrt(TTM);

% Compute Black's digital price 
bl_price=DigitalPriceBlack(F0, K, digital_vol, TTM, expirydiscount, notional);

%compute the price considering the digital risk
digital_price = bl_price-slope_imp.*vega.*notional;

end