% Author: Zi Wang
function [f, f2]= nnet_boston(nnsize, mu, mu_dec, mu_inc)
% The goal is to maximize the output f. (f = 20 - validation error)
% nnsize \in [1 100] 
% mu \in [0.000001, 100]
% mu_dec \in [0.01, 1]
% mu_inc \in [1.01, 20]
s = RandStream('mcg16807','Seed', 1024);
RandStream.setGlobalStream(s);

load boston_data.mat
nnsize = floor(nnsize);
net = feedforwardnet(nnsize);
net.trainParam.max_fail = 10;
net.trainParam.min_grad = 10^(-15);
net.trainParam.mu = mu;
net.trainParam.mu_dec = mu_dec;
net.trainParam.mu_inc = mu_inc;
net.divideFcn = 'divideint';
net.divideParam.trainRatio = 0.8;
net.divideParam.valRatio = 0.1;
net.divideParam.testRatio = 0.1;
net.trainParam.time = 120;
[net,tr] = train(net,x,t);
y = net(x);
f2 = perform(net,y,t);
f = 20-tr.best_vperf; 