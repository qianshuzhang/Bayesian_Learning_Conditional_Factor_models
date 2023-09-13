function X = PriorSim(param, Nsample, n_factors)

% Parameters need to be estimated:
% 1.mu, 2.omega, 3.alpha, 4.beta

X = param.mu + param.sig.*randn(Nsample,1);
for i = 1:(n_factors*4)
    if i>n_factors
        prior = param.mu(i)+param.sig(i).*randn(Nsample*10,1);
        j = find(prior > 0);
        j = j(randperm(length(j)));
        X(:, i) = prior(j(1:Nsample));
    else
        X(:,i) = param.mu(i)+param.sig(i).*randn(Nsample,1);
    end
end

for i = (n_factors*4+1):(n_factors*4+(n_factors-1)*n_factors/2)
    prior = param.mu(i)+param.sig(i).*randn(Nsample*10,1);
    j = find(prior < 1 & prior >0);
    j = j(randperm(length(j)));
    X(:, i) = prior(j(1:Nsample));
end

% X(:,1:n_factors) = prior(1:Nsample,1:2);
% X = prior(1:Nsample,:);
% for i = (n_factors+1):(4*n_factors)
%     j = find(prior(:,i)>0);
%     X(:,i) = prior(j(1:Nsample),i);
% end






% X(:, 1) = param.mu(1) + param.sig(1)*randn(Nsample, 1);
% 
% prior = param.mu(2) + param.sig(2)*randn(Nsample*10, 1);
% j = find(prior > 0);
% X(:, 2) = prior(j(1:Nsample));
% 
% prior = param.mu(3) + param.sig(3)*randn(Nsample*10, 1);
% j = find(prior > 0);
% X(:, 3) = prior(j(1:Nsample));
% 
% prior = param.mu(4) + param.sig(4)*randn(Nsample*10, 1);
% j = find(prior > 0);
% X(:, 4) = prior(j(1:Nsample));





%{
    j = find(prior(:,n_factors*i+1)>0);
    X(:,n_factors*i+1) = prior(j(1:Nsample),n_factors*i+1);
    j = find(prior(:,n_factors*i+2)>0);
    X(:,n_factors*i+2) = prior(j(1:Nsample),n_factors*i+2);
%} 