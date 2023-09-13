function [l, lnw, t, States, Statistics] = ProcessNewobs(X, l, lnw, t, States, Statistics, Ret, Tpred,n_factors)

T = length(Ret);

Nparam = size(X, 1);

fnames_states = fieldnames(States);

ESS = Nparam;

while (t <= T & ESS > 0.5*Nparam)
    
    [incrementalw, States] = filtering_llh(Ret(t,:), X, States,n_factors);
    
    j = isnan(incrementalw) | imag(incrementalw) ~= 0;
        
    incrementalw(j) = -inf;
    
    l = l + incrementalw;
    
    % get incremental normalizing ratio
    W_prev = exp(lnw - max(lnw));
    W_prev = W_prev/sum(W_prev);
    
    max_incrementalw = max(incrementalw);
    
    normfac = log(sum(W_prev.*exp(incrementalw - max_incrementalw))) + max_incrementalw;
    
    Statistics.normfac = [Statistics.normfac; normfac];
    
    lnw = lnw + incrementalw;
    
    % normalized weights and ESS
    W = exp(lnw - max(lnw));
    W = W/sum(W);
    
    ESS = 1/sum(W.^2);
    
    Statistics.ESS = [Statistics.ESS; ESS];
    
    % volatility forecasting
    hfm = zeros(n_factors, Tpred);
    
    sigmaR2 = X(:, (n_factors*1+1):(n_factors*2))./(1 - X(:, (n_factors*2+1):(n_factors*3)) - X(:, (n_factors*3+1):(n_factors*4)));
    
    vf1 = X(:, (n_factors*1+1):(n_factors*2)) + X(:, (n_factors*2+1):(n_factors*3)).*States.e.^2 + X(:, (n_factors*3+1):(n_factors*4)).*States.h;
    
    hfm(:,1) = mean(vf1);
    
    for tn = 2:Tpred
        
        vfn = sigmaR2 + (X(:, (n_factors*2+1):(n_factors*3)) + X(:, (n_factors*3+1):(n_factors*4))).^(tn-1).*(vf1 - sigmaR2);
        
        hfm(:,tn) = mean(vfn);
        
    end
    
    Statistics.EVm = [Statistics.EVm; hfm];
    
    % States statistics
    for i = 1:length(fnames_states)
       
        DSeven = Resample_vec(States.(fnames_states{i}), W, 10*Nparam);
        
        DS_m(:,i)     = mean(DSeven);
        DS_prc5(:,i)  = prctile(DSeven, 5);
        DS_prc95(:,i) = prctile(DSeven, 95);
        
    end
    
    Statistics.DS_m     = [Statistics.DS_m; DS_m];
    Statistics.DS_prc5  = [Statistics.DS_prc5; DS_prc5];
    Statistics.DS_prc95 = [Statistics.DS_prc95; DS_prc95];
    
    % statistics for parameters
    Xeven = Resample_vec(X, W, 10*Nparam);
    
    Statistics.Xm      = [Statistics.Xm; mean(Xeven)];
    Statistics.X_prc5  = [Statistics.X_prc5; prctile(Xeven, 5)];
    Statistics.X_prc95 = [Statistics.X_prc95; prctile(Xeven, 95)];
    
    t = t + 1;
    
end