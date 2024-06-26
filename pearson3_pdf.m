function pdf = pearson3_pdf(x,alpha,beta,xi,rho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the probability density function for the Pearson 3 distribution 
% at the values in x. The Hosking and Wallis (1997) version of the 
% distribution is chosen.  
%
% Given the location (mu), scale (sigma) and shape (rho) parameters, the 
% three other parameters used in this version are
%   alpha = 4 / sigma^2
%   beta = 0.5 * sigma * abs(rho)
%   xi = mu - 2*sigma/rho
%
% If rho > 0, the range of x is : xi <= x < Inf 
% If rho = 0, the range of x is : -Inf < x < Inf 
% If rho < 0, the range of x is : -Inf < x <= xi 
%
% Input arguments
%    x                      vector of values
%    alpha, beta, xi        parameters of the distribution
%    rho                    shape parameter. Depending on its value, the 
%                           skewness is positive (rho > 0) or negative 
%                           (rho < 0). If rho = 0, the distribution is 
%                           normal where the mean is alpha and the standard
%                           deviation is beta
% Output argument(s)
%   pdf                     vector of probability densities of PIII(xi,a,b)
%
% Reference: 
%   Hosking, J., & Wallis, J. (1997). Regional Frequency Analysis:
%       An Approach Based on L-Moments. Cambridge: Cambridge University Press. 
%       doi:10.1017/CBO9780511529443
%
% Guillaume Talbot, INRS-ETE 2021
% modified by Jasper A. Vrugt, July 2021
% UCI
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if rho == 0 % No skewness, treat as N(mu,s2)
    x2 = (x-alpha)./beta;
    pdf = exp(-0.5.*x2.^2)./sqrt(2.*pi)./beta;
elseif rho > 0 % Positive skewness    
    pdf = (x-xi).^(alpha-1).*exp((xi-x)./beta)./beta.^alpha./gamma(alpha);
elseif rho < 0 % Negative skewness 
    pdf = (xi-x).^(alpha-1).*exp((x-xi)./beta)./beta.^alpha./gamma(alpha);
end

% When the pdf is complex (i.e. out of range), replace it with zeros
pdf(imag(pdf)~=0) = 0;
% JAV: I do not like this, but OK for now [inconsequential in our use]