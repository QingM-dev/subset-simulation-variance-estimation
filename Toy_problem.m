%% =======================================================
%  Toy example: variance of correlated Gaussian sequences
%  and proposed variance estimation by 'moving block'
% =======================================================

clear; clc;

%% --- Settings ---
nSet   = 1e4;      % number of realizations
nChain = 1000;     % number of chains
rhoVal = 0.7;      % correlation decay parameter

rng('shuffle');

%% --- Build correlation matrix (Toeplitz) ---
lags = rhoVal .^ (0:nChain-1);
C    = toeplitz(lags);           % correlation matrix
mu   = zeros(nChain, 1);

%% --- Simulate correlated chains ---
X    = zeros(nSet, nChain);
sVar = zeros(nSet, 1);           % sample variance for each chain

for i = 1:nSet
    X(i,:) = mvnrnd(mu, C);
    sVar(i) = var(X(i,:));
end

%% --- Insample variance estimation by 'moving block' method ---
nMax  = nChain/4;      % max block size limit (can be adjusted)
nUse  = 5000;          % number of realizations
Rmat  = zeros(nUse, nMax);
vEst  = zeros(nUse, 1);

for i = 1:nUse
    
    vBlk      = zeros(nMax,1);
    vBlk(1)   = sVar(i);
    
    for b = 2:nMax
        % block means
        B = movmean(X(i,:), b, 'Endpoints', 'discard');
        
        % variance of block means
        vBlk(b) = var(B);
        
        % R factor
        Rmat(i,b) = b * vBlk(b) / vBlk(1);
    end
    
    % Reff  (exclude b=1)
    Reff = max(Rmat(i,2:end));
    
    % variance estimate of sample mean
    vEst(i) = sVar(i) * Reff / nChain;
end

%% --- “Actual” variance based on all data ---
mSet = mean(X,2);               % mean of each chain
vAll = var(X, 0, 'all');        % variance of all entries

Rtrue = var(mSet) * nChain / vAll;

Rmax_all = max(Rmat(:,2:end), [], 2);


%% --- Comparison plots: estimated vs “actual” variance of mean ---
figure;
histogram(vEst, 40, 'Normalization', 'pdf'); hold on;
xline(var(mSet), 'LineWidth', 2);
hold off;

xlabel('Variance of sample mean');
ylabel('PDF (estimated)');
title('Estimated Var(\bar{X}) from block method vs empirical truth');
legend('Estimated Var(\bar{X})', 'Empirical Var(\bar{X})', 'Location', 'best');
grid on;

