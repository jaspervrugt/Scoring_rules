function X = pearson3_rnd(alpha,beta,xi,rho,m,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute a random sample of m x n values from the Pearson 3 distribution. 
% The Hosking and Wallis (1997) version of the distribution is chosen.  
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
%    alpha, beta, xi        parameters of the distribution
%    rho                    shape parameter. Depending on its value, the 
%                           skewness is positive (rho > 0) or negative 
%                           (rho < 0). If rho = 0, the distribution is 
%                           normal where the mean is alpha and the standard
%                           deviation is beta
%    m, n                   number of rows and columns of matrix X 
% Output argument(s)
%    X                      m x n matrix of random samples from PIII(xi,a,b)
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

% Generate a m x n matrix with values between 0 and 1
p = rand(m,n);

% Draw a matrix X of size m x n with samples from PIII(xi,a,b)
if rho == 0 % No skewness, draw from normal distribution
    X = norminv(p,alpha,beta);    
elseif rho > 0 % Positive skewness    
    X = gammaincinv(p,alpha.*ones(size(p)));
    X = X.*beta + xi;  
else % Negative skewness    
    X = gammaincinv(1-p,alpha.*ones(size(p)));    
    X = - X.*beta + xi;
end