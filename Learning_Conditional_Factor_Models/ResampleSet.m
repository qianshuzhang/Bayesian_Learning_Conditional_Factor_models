function [X,lnw,l,Particles,lnpdf_prior]=ResampleSet(X,lnw,l,Particles,lnpdf_prior,n_factors)
%this function resamples the fixed parameter set


Nparam=size(X,1);

cumw = [0; cumsum(exp(lnw - max(lnw)))];
cumw = cumw/cumw(end);

%use stratified resampling
[PH, bin] = histc(((0:Nparam-1)+rand(1,Nparam))/Nparam,cumw);
X = X(bin,:);
l = l(bin);

lnpdf_prior = lnpdf_prior(bin);

FNames  = fieldnames(Particles);

Nfields = length(FNames);

for i = 1:Nfields
    Particles.(FNames{i}) = Particles.(FNames{i})(bin,:);
end

lnw = zeros(Nparam, 1);


%{
1. Nparam=size(X,1);：这行代码确定了粒子或元素的数量，其中X是一个数据矩阵，size(X,1)返回矩阵的行数。

2. cumw = [0; cumsum(exp(lnw - max(lnw)))];：
这行代码计算了经过归一化处理的累积权重。lnw是一个包含粒子权重的向量，
exp(lnw - max(lnw))首先对权重进行归一化处理，然后计算指数，cumsum函数计算累积和。
[0; cumsum(exp(lnw - max(lnw)))]将0添加到累积和的开头，确保第一个区间的边界。

3. cumw = cumw/cumw(end);：这行代码将累积权重归一化，使其总和等于1。
通过将累积权重除以最后一个元素，确保归一化后的累积和是1。

4. histc(((0:Nparam-1)+rand(1,Nparam))/Nparam, cumw);：这行代码执行了分层重采样的关键步骤。
首先，生成一个介于0和1之间的随机数向量，(0:Nparam-1)+rand(1,Nparam)生成了一个包含Nparam个随机数的向量，范围在[0, 1)之间。
然后，将这个向量除以Nparam，将其范围缩放到[0, 1/Nparam)之间。
最后，使用histc函数根据归一化的累积权重cumw，将这些随机数映射到对应的区间，并返回每个随机数所属的区间索引。
这样就完成了分层重采样，根据权重进行抽样并获得对应的区间索引。

5. X = X(bin, :);和l = l(bin);：这两行代码根据分层重采样得到的区间索引bin，
从数据矩阵X和标签向量l中选择对应的样本。这样，通过分层重采样，根据权重重新选择了样本，
以获得更平衡和代表性的样本集合。
%}