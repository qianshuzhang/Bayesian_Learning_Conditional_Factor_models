function [X_proposal, lnprop_X, lnprop_Xprop] = SimProposalTrMixtureNormal(X, Ndim, Ntr,n_factors,vary_param)

%simulate from mv normal fitted on X allows transformed marginals
%(exp,logistic)
%assume even-weighted sample

Nparam = size(X, 1);

X_vary = X(:,vary_param);

mixturedim = Ndim;

scaleX = repmat(std(X_vary), Nparam, 1);

X_vary = X_vary./scaleX;

X_tr = X_vary;

X_tr_proposal = X(1,:).*ones(Ntr*Nparam,1);

% fit mv mixture normal on X_tr and simulate from it

options = statset('Display','off','MaxIter', 500);

try
    
    obj = gmdistribution.fit(X_tr, mixturedim, 'Options', options);
    
catch
    
    obj = gmdistribution.fit(X_tr, mixturedim, 'SharedCov', true, 'Options', options);
    
end

X_tr_proposal(:,vary_param) = random(obj, Ntr*Nparam);

lnprop_Xprop = log(pdf(obj, X_tr_proposal(:,vary_param)));

lnprop_X = log(pdf(obj, X_tr));

% transform proposal

X_proposal = X_tr_proposal;

scaleX = repmat(scaleX, Ntr, 1);

X_proposal(:,vary_param) = X_proposal(:,vary_param).*scaleX;

