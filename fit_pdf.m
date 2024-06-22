function sr_CV = fit_pdf(x,l,CV_opt,gPsetting)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Function to help fit the specified distributions to given CV                       %%
%%                                                                                    %%
%% SYNOPSIS: sr_CV = fit_pdf(x,opt,CV_opt,gPsetting)                                  %%
%%  where                                                                             %%
%%   x         [input]  REQUIRED: vector of shape parameters distribution             %%
%%   l         [input]  REQUIRED: integer of distribution to use/fit                  %%
%%   CV_opt    [input]  REQUIRED: desired value of coefficient of variation           %%
%%   gPsetting [input]  REQUIRED: value scale parameter (= sigma) of GEV distribution %%
%%   sr_CV     [output] squared residual of CV                                        %%
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt                                                     %%
%% University of California Irvine                                                    %%
%% July 2021                                                                          %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %% 

% Lazy approach: we draw N samples from specified distribution and then use 
% the sample mean and sample variance to determine the CV
% Better approach: Compute mean and variance directly from values of shape
% parameters. This is more accurate but inconsequential

N = 1e5; % Number of samples to draw from specified distribution
switch l
    case 2 % normal distribution
        X = normrnd(x(1),x(2),N,1);
    case 3 % Generalized extreme value distribution
        X = gevrnd(x(1),1,x(2),N,1);
    case 4 % lognormal distribution
        X = lognrnd(x(1),x(2),N,1);
    case 5 % Generalized Pareto distribution
        X = gprnd(x(1),gPsetting,x(2),N,1); 
end
mX = mean(X);                   % Sample mean       
sX = std(X);                    % Sample standard deviation    
CV = sX/mX;                     % Sample CV    
sr_CV = sum((CV_opt - CV).^2);  % squared CV residual

end