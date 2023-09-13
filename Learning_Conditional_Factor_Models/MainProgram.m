% Main program: Learning Comovements of Tails, Volatility and Returns

clc;
clear;

%% Data simulation or data loading
load('SP500_Jan80_Dec12');
Ret = diff(log(SP500_Jan80_Dec12(:, 2)));
T   = length(Ret);
    
%% ------------implement Online Learning-----------------------------------
% -------------------------------------------------------------------------
%particle numbers
Nparam    = 2*1024;
ESS_Bound = Nparam/2;

% Parameters need to be estimated:
% [1.mu, 2.omega, 3.alpha, 4.beta];
% among which omega, alpha, beta need to be positive

prior_param.mu  = [0   , 1e-6, 0.05, 0.90];
    
prior_param.sig = [2e-4, 2e-6, 0.05, 1.20];

% simulate from the prior of the fixed params
X = PriorSim(prior_param, Nparam);

lnpdf_prior  = PriorLogl(X, prior_param);

l   = zeros(Nparam, 1);
lnw = zeros(Nparam, 1);

t = 1; 

%initial values of volatilities
States.h  = 0.025/252 * ones(Nparam, 1);
States.e  = sqrt(States.h).*randn(Nparam, 1);

%initialize statistics structure
Statistics.ESS = [];
Statistics.normfac = [];

% will store posterior mean of dynamic quantities
Statistics.DS_m     = [];
Statistics.DS_prc5  = [];
Statistics.DS_prc95 = [];

% will store posterior mean and 5,95 percentiles for fixed parameters
Statistics.Xm      = [];
Statistics.X_prc5  = [];
Statistics.X_prc95 = [];

% volatility forecasting
Statistics.EVm = [];

% counter of move steps
movecounter = 0;

% structure to store runtime of processing steps
Runtime.Process_t = [1];
Runtime.Process_runtime = [0];

starttime = tic;
Runtime.total_runtime = toc(starttime);

Tpred = 7;

% begin to loop over the data
%%
while (t <= T)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % process new observations %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tic
    
    [l, lnw, t, States, Statistics] = ProcessNewobs(X, l, lnw, t, States, Statistics, Ret, Tpred);
    
    proctime = toc;
    
    Runtime.Process_t = [Runtime.Process_t, t];
    
    Runtime.Process_runtime = [Runtime.Process_runtime, proctime];
    
    avg_time = Runtime.Process_runtime(end)/(Runtime.Process_t(end) - Runtime.Process_t(end-1));
    
    disp(' ');
    
    disp(['Observations ' num2str(Runtime.Process_t(end-1)) '-' num2str(Runtime.Process_t(end))...
          ' are processed, at per observation time of ' num2str(avg_time) ' sec' ]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % resample and move if necessary %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if Statistics.ESS(end) < 0.5*Nparam
        
        disp(' ');
        disp(['Resample and move at ' num2str(t)]);
        
        %%%%%%%%%%%%
        % Resample %
        %%%%%%%%%%%%
        
        [X, lnw, l, States, lnpdf_prior] = ResampleSet(X, lnw, l, States, lnpdf_prior);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % while unique particles are less than 0.5*Nparam, move %
        % using a mixture normal independent proposal                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        counter = 0;
        
        while or(counter == 0, size(unique(X, 'rows'), 1) < 0.5*Nparam)
            
            counter = counter + 1;
            
            movecounter = movecounter + 1; Moves.t(movecounter) = t;
            
            tic
            
            [X, States, l, lnpdf_prior, Moves.AcceptRate(movecounter)] = ...
                    MoveSet(X, States, l, lnpdf_prior, prior_param, Ret(1:t));
            
            Moves.Runtime(movecounter) = toc;
            
            disp(' ');
            
            disp( ['Move number ' num2str(counter) ' at t=' num2str(t) ' took ' num2str(Moves.Runtime(movecounter)) ' sec']);
            
            disp(['unique X after move number ' num2str(counter) ' : ' num2str(size(unique(X,'rows'),1))]); 
        end
    end
    
    Runtime.total_runtime = [Runtime.total_runtime, toc(starttime)];
    
    save estresults_GARCH Statistics t X lnw Moves Runtime;
    
end
