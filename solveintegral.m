function I=solveintegral(x, integrand,extraparam)
% solventegral solves the fourier transform of a given integrand function.

% I = solveintegral(x, integrand) computes the integral function of the given 
% integrand function using numerical integration and evaluates it over the range    
% specified by the vector x .

% I = solveintegral(x, integrand, extraparam) computes the integral
% using FFT integration with the provided extrapolation parameters and valuates 
% it over the range  specified by the vector x.

% Inputs:
% x:            Vector of evaluation points for the integral.
% integrand:    Function handle representing the integrand.
% extraparam:   Structure containing FFT initialization parameters.
%               - For initialization with dz: extraparam contains
%                 fields M and dz.
%               - For initialization with psi_start: extraparam
%                 contains fields M and psi_start.

% Outputs:
% I: Vector of computed integrals corresponding to the evaluation
%    points in x.

% Initialize integral vector
I=zeros(1, length(x));

switch nargin()
    
    % FFT algorithm
    case 3
        %Initialization
        N=2^(extraparam.M);
        
        % Case where I initialized dz
        if isempty(extraparam.psi_start)
            dz=extraparam.dz;
            z_start = -(N-1)*dz/2;
            dpsi=2*pi/(N*dz); 
            psi_start = -(N-1)*dpsi/2;
        
        % Case where I initialized psi_start
        else if isempty(extraparam.dz)
            psi_start=extraparam.psi_start; 
            dpsi = -2*psi_start/(N-1);  
            dz = 2*pi/(N*dpsi);  
            z_start = -(N-1)*dz/2;  
            psi_start = -(N-1)*dpsi/2;
        
        % Return error if neither are initialized
        else
            error("Non valid initialization method")
        end
        end

        % Compute the integral using FFT algorithm

        % Row vector of integration points
        psi=linspace(psi_start, -psi_start, N);

        % Row vectors containing moneyness evaluation points
        z=linspace(z_start, -z_start, N);

        % Row vector of rotated integrand valuated on integration points 
        rotintegrand=integrand(psi).*exp(-1i.*[0:N-1]*dpsi*z_start);
        
        % Row vector of prefactors
        prefactor=dpsi*exp(-1i*psi_start.*z);

        % Transformed integrand valuated at moneyness points 
        I=real(prefactor.*fft(rotintegrand));

        % Interpolate the integral through linear interpolation to evaluate
        % it in the give points
        I=interp1(z,I,x,"linear");

    % Use built-in MATLAB quadrature methods
    case 2
        
        % For each x define the integrand and compute the integral
        for ii=1:length(x)
            tot_integrand = @ (psi) exp(-1i* x(ii)* psi).* integrand(psi);
            I(ii) = real(quadgk(tot_integrand, -inf, inf));
        end

    % If nor 2 or 3 parameters were passed to the function: error
    otherwise
        error("Integration method not available")

end
end