function d_KL = KL_div(P,Q,b,KLoption)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Küllback-Leibler divergence of two discrete probability distributions, P and Q     %%
%%                                                                                    %%
%% SYNOPSIS: d_KL = KL_div(P,Q,b,KLoption)                                            %%
%%  where                                                                             %%
%%   P         [input]  REQUIRED: n x m matrix of probabilities (m different P's)     %%
%%   Q         [input]  REQUIRED: n x m matrix of probabilities (m different Q's)     %%
%%   b         [input]  OPTIONAL: base of logarithm [default = 2]                     %%
%%   KLoption  [input]  OPTIONAL: [0] Forward and [1] Forward+reverse KL divergence   %%
%%   d_KL      [output] 1 x m vector with forward KL-divergence, d(P||Q) or           %%
%%                      2 x m vector with forward d(P||Q) and reverse KL_div, d(Q||P) %%
%%                                                                                    %%
%% Based on download of KLDiv from Mathworks                                          %%
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt                                                     %%
%% University of California Irvine                                                    %%
%% July 2021                                                                          %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Check number of input arguments
if nargin < 4
    KLoption = 0;   % Only return forward KL divergence, D(P||Q)
end
if nargin < 3
    b = 2;          % Default base of the logarithm --> unit of bits
end
[n,m] = size(P);
% Check dimension of input arguments
if size(Q,1) ~= n
    error('KL_div: The number of rows of Q should match number of rows of P');
end
if size(Q,2) ~= m
    error('KL_div: The number of columns of Q does not match counterpart of P');
end
% Check values of input arguments
if sum(isinf(P))
    error('KL_div: At least one entry of matrix P is infinite!')
end
if sum(isinf(Q))
    error('KL_div: At least one entry of matrix Q is infinite!')
end
if any(P < 0)
    error('KL_div: At least one entry of matrix P is negative!')
end
if any(Q < 0)
    error('KL_div: At least one entry of matrix Q is negative!')
end
% Normalize P and Q
p = P./sum(P); q = Q./sum(Q);
% Compute product, p(x) * log(p(x)/q(x))
E = p.*log2(p./q); 
% Limit: p*log(p) = 0 if p goes to 0.
idx = (p < realmin); E(idx) = 0; 
% Return forward KL divergence, D(P||Q)
d_KL = sum(E,1);
% Check whether to compute reverse divergence, D(Q||P)
if KLoption == 1
    % Compute product, q(x)*log(q(x)/p(x))
    E = q.*log2(q./p);
    % Limit: q*log(q) = 0 if q goes to 0.
    idx = (q < realmin); E(idx) = 0; 
    % Return forward KL divergence, D(P||Q)
    d_KL(2,1:m) = sum(E,1); 
end
% Now adjust to specified base b
d_KL = d_KL/log2(b);

end