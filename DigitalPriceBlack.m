function digital_price=DigitalPriceBlack(F0, K, volatility, TTM, expirydiscount, notional)
% DigitalPriceBlack: computes Black's price of a digital option

% Inputs:
% F0:               initial forward
% K:                strike
% volatility:       volatility of the underlyine
% TTM:              time to maturity of the option
% expirydiscount:   discount at expiry (considering TTM)
% notional:         notional 

% Outputs:
% digital_price:    price of the digital option

% Compute d2 of the black model
d2=1./sqrt(TTM.*volatility.^2).*(log(F0./K) - 1/2.*TTM.*volatility.^2);

% Apply Black's formula
digital_price=expirydiscount.*normcdf(d2).*notional;

end