function lnpdf_prior = PriorLogl(X, prior_param,n_factors)

% 1.mu, 2.omega, 3.alpha, 4.beta

Nsample = size(X, 1);


E = repmat(prior_param.mu, Nsample, 1);
V = repmat(prior_param.sig, Nsample, 1);

lnpdf_prior = sum(-(X - E).^2./(2*V.^2) - log(2*pi)/2 - log(V), 2);


% enforce parameters in their supports; 2, 3, 4
for i = 1:n_factors
    j = X(:, i+1*n_factors) < 0 | X(:, i+2*n_factors) < 0 | X(:,i+3*n_factors) < 0 | X(:, i+2*n_factors) + X(:, i+3*n_factors) >= 3;
    lnpdf_prior(j) = -inf;
end

for i = (n_factors*4+1):(n_factors*4+(n_factors-1)*n_factors/2)
    j = X(:,i)<0 | X(:,i) >1;
    lnpdf_prior(j) = -inf;
end

% lnpdf_prior = sum(lnpdf_prior,2);