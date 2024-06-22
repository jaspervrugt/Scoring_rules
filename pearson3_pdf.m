%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the probability density function for the Pearson 3
% distribution at each of the values in x. The Hosking and Wallis (1997) 
% version of the distribution is chosen.  
%
% Given the location (mu), scale (sigma) and shape (Gamma) parameters, we can
% estimate the three other parameters used in this version :
%   -alpha : 4 / sigma^2
%   -beta :0.5 * sigma * abs(Gamma)
%   -xi : mu - 2*sigma/Gamma
%
% If Gamma  > 0, the range of x is : xi <= x < Inf 
% If Gamma  = 0, the range of x is : -Inf < x < Inf 
% If Gamma  < 0, the range of x is : -Inf < x <= xi 
%
% Input :
%    -x : vector of values
%    -alpha, beta, xi : parameters of the distribution
%    -Gamma : shape parameter. Depending on its value, the skewness is
%       positive (Gamma > 0) or negative (Gamma < 0). If Gamma = 0, the
%       distribution is normal where the mean is alpha and the standard
%       deviation is beta
%
% Output
%   -pdf : vector of the probability density function
%
% Source : Hosking, J., & Wallis, J. (1997). Regional Frequency Analysis:
% An Approach Based on L-Moments. Cambridge: Cambridge University Press. 
% doi:10.1017/CBO9780511529443
%
% Guillaume Talbot, INRS-ETE 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [pdf]=pearson3_pdf(x,alpha,beta,xi,Gamma)

if Gamma>0 %Case of positive skewness    
    pdf=(x-xi).^(alpha-1).*exp((xi-x)./beta)./beta.^alpha./gamma(alpha);
elseif Gamma<0 %Case of negative skewness 
    pdf=(xi-x).^(alpha-1).*exp((x-xi)./beta)./beta.^alpha./gamma(alpha);
else %Case Normal distribution (i.e Gamma = 0)
    x2=(x-alpha)./beta;
    pdf=exp(-0.5.*x2.^2)./sqrt(2.*pi)./beta;
end

%When the pdf is complex (i.e. out of range), replace it with zeros
pdf(imag(pdf)~=0)=0;