function [X, States, l, lnpdf_prior, AcceptRate] = MoveSet(X, States, l, lnpdf_prior, prior_param, Ret,initial_h,n_factors,vary_param)


Nparam = size(States.h, 1);

fnames_states = fieldnames(States);

T = length(Ret);

% set proposal density

Ndim = 3;
Ntr  = 20;


[X_proposal, logprop_orig, logprop_proposal] = SimProposalTrMixtureNormal(X, Ndim, Ntr,n_factors,vary_param);
    
% compute prior at proposal
lnpdf_prior_proposal = PriorLogl(X_proposal, prior_param,n_factors);

% only run the filter for sets when prior is finite
keep = find(isfinite(lnpdf_prior_proposal));

lnpdf_prior_proposal = lnpdf_prior_proposal(keep(1:Nparam));

X_proposal = X_proposal(keep(1:Nparam), :);

logprop_proposal = logprop_proposal(keep(1:Nparam));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimated loglikelihood at proposal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

individual_l_proposal = -inf*ones(Nparam, T);

% n factors
States_proposal.h  = initial_h;
States_proposal.hf = 0.025/252 * ones(Nparam, 1);
States_proposal.e  = sqrt(States_proposal.h).*randn(Nparam, n_factors);

T = size(Ret,1);

for t = 1:T
    
    [individual_l_proposal(:,t), States_proposal] = filtering_llh(Ret(t,:), X_proposal, States_proposal,n_factors);
    
end

l_proposal  = sum(individual_l_proposal, 2);

% M-H steps

logtargetn_proposal = l_proposal + lnpdf_prior_proposal;

drop = imag(logtargetn_proposal) ~= 0;

logtargetn_proposal(drop) = -inf;

logtargetn_orig = l + lnpdf_prior;

lnalpha = logtargetn_proposal - logtargetn_orig + logprop_orig - logprop_proposal;

logu   = log(rand(Nparam, 1));

accept = find(logu < lnalpha);

AcceptRate = length(accept)/Nparam;

disp(['Acceptance rate in MH step: ' num2str(AcceptRate)]);

% implement move

X(accept, :) = X_proposal(accept, :);

l(accept) = l_proposal(accept);

lnpdf_prior(accept) = lnpdf_prior_proposal(accept);

for i = 1:length(fnames_states)
    
    States.(fnames_states{i})(accept,:) = States_proposal.(fnames_states{i})(accept,:);
    
end