function [X_proposal, lnprop_X, lnprop_Xprop] = SimProposalTrMixtureNormal(X, Ndim, Ntr,n_factors)

%simulate from mv normal fitted on X allows transformed marginals
%(exp,logistic)
%assume even-weighted sample

Nparam = size(X, 1);

mixturedim = Ndim;

scaleX = repmat(std(X), Nparam, 1);

X = X./scaleX;

X_tr = X;

% fit mv mixture normal on X_tr and simulate from it

options = statset('Display','off','MaxIter', 500);

try
    
    obj = gmdistribution.fit(X_tr, mixturedim, 'Options', options);
    
catch
    
    obj = gmdistribution.fit(X_tr, mixturedim, 'SharedCov', true, 'Options', options);
    
end

X_tr_proposal = random(obj, Ntr*Nparam);

lnprop_Xprop = log(pdf(obj, X_tr_proposal));

lnprop_X = log(pdf(obj, X_tr));

% transform proposal

X_proposal = X_tr_proposal;

scaleX = repmat(scaleX, Ntr, 1);

X_proposal = X_proposal.*scaleX;

