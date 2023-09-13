function [l, States] = filtering_llh(Rt, X, States,n_factors)

% INPUTS:
%    X    : Nparam x K, where K is the number of parameters, and Nparam is
%            the number of parameter sets (param particles)
%            Nparam = size(params, 1);     % the number of parameter sets
%         
%    state: .e, .h1, Nparam x 1
%
% OUTPUTS:
%    l : Nparam x 1 (likelihood for each parameter set at time t)
% state: .e, .h, .hf   Nparam x 1

params.mu    = X(:, 1:n_factors);
params.omega = X(:, (n_factors*1+1):(n_factors*2));
params.alpha = X(:, (n_factors*2+1):(n_factors*3));
params.beta  = X(:, (n_factors*3+1):(n_factors*4));

lower_tri_elements =  X(:, (n_factors*4+1):end);
for dim = 1:size(X,1) 
    M = zeros(n_factors); % n维
    
    % 使用逻辑索引将下三角元素填充回去
    M(tril(true(n_factors), -1)) = lower_tri_elements(dim,:);
    
    % 将对角线元素设置为1
    M = M + eye(n_factors);
    
    params.A(:,:,dim) = M;
end

e = States.e;
h = States.h;

h  = max(params.omega + params.beta.*h + params.alpha.*e.^2, 1e-25);
  
mu = params.mu;

pdf_values = zeros(n_factors,1);
for i = 1:size(X,1)
    sigma = params.A(:,:,i)*diag(h(i,:))*(params.A(:,:,i)');
    pdf_values(i) = -log(2*pi) *(n_factors/2) - log(det(sigma))/2- 0.5 * (Rt - mu(i,: )) * inv(sigma) * (Rt - mu(i, :))';
end

l = pdf_values;
% l = sum(-(Rt - mu).^2./(2*h) - log(h)/2 - log(2*pi)/2,2);

e = Rt - mu;

States.e  = e; 
States.h  = h;


