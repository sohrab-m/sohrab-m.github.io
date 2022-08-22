%% Prediction with Noise-free Observations
%
% Generates Figure 3. Contains random functions, expect values to vary. If
% reproducibility is required - specify some random seed.
%

%%
clear all
close all
colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k'];
n_sequences = 3; % number of sequences/functions 
%-----------------------------------
x_lim = [-5, 5];

n = 50;
X_true = linspace(x_lim(1), x_lim(2), n);
f_true = sin(X_true); % model is a sin funciton

mu = zeros(n, 1);
K_star_star = cov_se(X_true, X_true);
sigma = K_star_star;
f_drawn = mvnrnd(mu, sigma, n_sequences);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphics (Prior)

figure(1)
stdev = sqrt(diag(K_star_star));
curve_p_stdev = mu + 2*stdev;
curve_m_stdev = mu - 2*stdev;
gtfill(X_true,curve_m_stdev,curve_p_stdev, 1000, 'm', [0, 0, 0]+0.9)
hold on
%plot(X_true, f_true, 'k--', 'linewidth', 0.5);
plot(X_true, f_drawn(1, :), '.', 'linewidth', 2, 'color', colors(1));

for seq = 2:n_sequences
    plot(X_true, f_drawn(seq, :), '-', 'linewidth', 2, 'color', colors(1+rem(seq-1, length(colors))));
end
xlabel('input, x');
ylabel('output, f(x)');
ylim([-5, 5]);


%% Posterior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_samples = 4; % num of noise-free observation points
X_sampled = sort(x_lim(1)+rand(1,n_samples)*(x_lim(2)-x_lim(1)));
f_sampled = sin(X_sampled);  % model is a sin funciton


K_xx = cov_se(X_sampled, X_sampled);
K_star_x = cov_se(X_true, X_sampled);
K_x_star = cov_se(X_sampled, X_true);
mu_new = K_star_x / K_xx * f_sampled';
sigma_new = K_star_star - K_star_x /K_xx * K_x_star;
%sigma_new = triu(sigma_new)+triu(sigma_new,1)';
f_new = mvnrnd(mu_new, sigma_new, n_sequences);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphics (Posterior)

figure(2)
stdev = sqrt(diag(sigma_new));
curve_p_stdev = mu_new + 2*stdev;
curve_m_stdev = mu_new - 2*stdev;
gtfill(X_true,curve_m_stdev,curve_p_stdev, 1000, 'm', [0, 0, 0]+0.9)

hold on
for seq = 1:n_sequences
    plot(X_true, f_new(seq, :), '-', 'linewidth', 2, 'color', colors(1+rem(seq-1, length(colors))));
end
xlabel('input, x');
ylabel('output, f(x)');
ylim([-5, 5]);
%plot(X_true, f_true, 'k--', 'linewidth', 0.5);
plot(X_sampled, f_sampled(1, :), 'k+', 'markersize', 20);


%% Extra
function K = cov_se(X, Y)
n = size(X, 2);
m = size(Y, 2);

K = zeros(n, m);
for i =1:n
    for j=1:m
        x = X(:, i);
        y = Y(:, j);
        K(i, j) = squared_exp(x, y);    
    end
end
end

function k = squared_exp(x, y)
k = exp(-1/2*norm(x-y)^2);
end
