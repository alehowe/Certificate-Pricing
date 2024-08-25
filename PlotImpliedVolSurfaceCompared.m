function PlotImpliedVolSurfaceCompared(F0, strikes, surface, TTM, expirydiscount, model_prices)
% PlotImpliedVolSurfaceCompared: plots the market and the model's implied
% volatilities surfaces given the model call prices

% Inputs:
% F0:                   initial forward
% strikes:              vector containing different strikes
% surface:              vector containing the volatility surface corresponding
%                       to the vector K
% TTM:                  time to maturity of the option
% expirydiscount:       discount at expiry (considering TTM)
% model_prices:         vector cotaining the model call prices at the given
%                       strikes defined in strikes above

% Compute the zero rate
r=-log(expirydiscount)/TTM;

% Initialize the calibrated implied volatility surface
implied_vol_model=zeros(1, length(strikes));

% Compute the calibrated implied volatility surface
implied_vol_model=blkimpv(F0, strikes, r, TTM, model_prices);

% Plot the market and the model's implied volatility surfaces
figure()
plot(strikes,surface, "-r");
hold on;
plot(strikes,implied_vol_model, "-b");
grid on;
legend("Market implied volatility surface", "Model implied volatility surface")
title("Implied volatility surface")
xlabel("Strikes")
ylabel("Impled volatilities")
