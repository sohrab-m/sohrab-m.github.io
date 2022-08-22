
clear;
clc;
close all;
%Get discretized dynamics
PlantModelQuadSimpleLinear;

rpmax = deg2rad(45); %max roll pitch angle
thrmax = 5*1; %max thrust

%Collect N trials (rollouts), horizon T, find A_hat, B_hat

N = 200;
T = 50;
x0= zeros(6,1);
initial_data_xt  = [];
initial_data_xt_1 = [];
std_u_rp = rpmax/3;
std_u_t = thrmax/3;
sigma_w = 0.1;
fail_prob = 0.01;

big_A = [sys_d.A zeros(6,3);zeros(6,6) sys_d.B];
E_A = zeros(N,1);
E_B = zeros(N,1);
bound_A_init = zeros(N,1);
bound_B_init = zeros(N,1);
bootstrap_eA = zeros(N,1);
bootstrap_eB = zeros(N,1);


for n_ind  = 1:N
    temp_x = zeros(T,9);
    temp_y = zeros(T,6);
    u0 = [std_u_rp*randn(2,1); std_u_t*randn(1,1)];
    z0= [x0 ; u0];
    x_t = x0;
    u_t = u0;
    z_t = z0;
    for t_ind = 1:T
        x_t_1 = sys_d.A*x_t +sys_d.B*u_t + sigma_w*randn(6,1);
        temp_y(t_ind, :) = x_t_1';
        temp_x(t_ind, :) = [x_t;u_t]';
        
        u_t = [std_u_rp*randn(2,1); std_u_t*randn(1,1)];
        x_t = x_t_1;
    end
    initial_data_xt = [initial_data_xt; temp_x];
    initial_data_xt_1 = [initial_data_xt_1; temp_y];
    
    
    beta = mvregress(initial_data_xt,initial_data_xt_1);
    A_hat  = beta(1:6, :)';
    B_hat  = beta(7:9, :)';
    
    
    %Compute initial error
    E_A(n_ind) = norm(A_hat - sys_d.A);
    E_B(n_ind) = norm(B_hat - sys_d.B);
    
    %Compute intial bounds
    C_squared = sigma_w^2*( sqrt(9)+ sqrt(6)+ sqrt(2*log(1/fail_prob)))^2;
    error= [(A_hat - sys_d.A)'; (B_hat - sys_d.B)'];
    error_norm = error*error';
    M_mat = inv(initial_data_xt' * initial_data_xt);
    error_bound_matrix = C_squared*M_mat;
    
    %Compute bound for A, B
    QA = [eye(6) zeros(6,3)];
    QB = [zeros(3,6) eye(3)];
    
    bound_A_init(n_ind) = sqrt(C_squared)*sqrt(norm(QA*M_mat*QA'));
    bound_B_init(n_ind) = sqrt(C_squared)*sqrt(norm(QB*M_mat*QB'));
    
    
    fprintf('Rollout %d\n',n_ind)
    fprintf('Initial E_A: %g\n', E_A(n_ind))
    fprintf('Initial Bound: %g\n', bound_A_init(n_ind))
    fprintf('Initial E_B: %g\n', E_B(n_ind))
    fprintf('Initial Bound: %g\n', bound_B_init(n_ind))
      
    %Boostrap estimation of E_A and E_B, M times
    
    M= 100;
    L = 50;
    bootstrap_data_xt = [];
    bootstrap_data_xt_1 = [];
    bootstrap_eA_set = zeros(M,1);
    bootstrap_eB_set = zeros(M,1);
    for m_ind  = 1:M
        for l_ind = 1:L
            temp_x = zeros(T,9);
            temp_y = zeros(T,6);
            u0 = [std_u_rp*randn(2,1); std_u_t*randn(1,1)];
            z0= [x0 ; u0];
            x_t = x0;
            u_t = u0;
            z_t = z0;
            for t_ind = 1:T
                x_t_1 = A_hat*x_t +B_hat*u_t + sigma_w*randn(6,1);
                temp_y(t_ind, :) = x_t_1';
                temp_x(t_ind, :) = [x_t;u_t]';
                
                u_t = [std_u_rp*randn(2,1); std_u_t*randn(1,1)];
                x_t = x_t_1;
            end
            bootstrap_data_xt = [bootstrap_data_xt; temp_x];
            bootstrap_data_xt_1 = [bootstrap_data_xt_1; temp_y];
        end
        
        %Find A_tilde, B_tilde
        
        beta = mvregress(bootstrap_data_xt,bootstrap_data_xt_1);
        A_tilde = beta(1:6, :)';
        B_tilde  = beta(7:9, :)';
        bootstrap_eA_set(m_ind) = norm(A_tilde - A_hat);
        bootstrap_eB_set(m_ind) = norm(B_tilde - B_hat);
    end
    
    %output 100(1-delta) percentile
    bootstrap_eA(n_ind) = quantile(bootstrap_eA_set, 1-fail_prob);
    bootstrap_eB(n_ind) = quantile(bootstrap_eB_set, 1-fail_prob);
    
    fprintf('Bootstrap E_A: %g\n', bootstrap_eA(n_ind))
    fprintf('Bootstrap E_B: %g\n', bootstrap_eB(n_ind))
    
end

save ('results.mat')

figure; plot(E_A); hold on; plot(bound_A_init); plot(bootstrap_eA);
legend({'E_A'; 'Initial Bound'; 'Bootstrap Bound'})
title('Error estimates vs iterations for A')
saveas(gcf, 'E_A_plot2.png')

figure; plot(E_B); hold on; plot(bound_B_init); plot(bootstrap_eB);
legend({'E_B'; 'Initial Bound'; 'Bootstrap Bound'})
title('Error estimates vs iterations for B')
saveas(gcf, 'E_B_plot2.png')



