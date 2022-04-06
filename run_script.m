%MC_EM for hawkes process param estimation with parameters v,alpha,beta 

%% Select if working univariate, bivariate or trivariate:
p = 2;
n_reps = 1; % how many times to simulate from the Hawkes ground truth



for case_idx = 1:n_reps %1:50 %%%13%49:500

    if p == 2
        disp('simulating data')
        %% Bivariate case
        % bin_width is the aggregation level. 
%         v = [1;1]; A = [0.9, 0; 0, 0.9]; B = [1.3, 1.3; 1.3, 1.3]; bin_width = 1; end_time = 2e3;  
        v = [1;1]; A = [0.9, 0; 0, 0.9]; B = [1.3, 1.3; 1.3, 1.3]; bin_width = 0.5; end_time = 1e3;  
%         v = [0.3; 0.3]; A = [0.7, 0.9; 0.6, 1]; B = [1.5, 2; 2, 3.5]; bin_width = 1; end_time = 2e3;

        % Check stationarity
        if p>1 && eigs(A./B, 1, 'lm')-1 >= 0 
            disp('non stationary parameters chosen')
            return
        elseif any(any(A>B)) && p == 1
            disp('non stationary parameters chosen')
            return
        end
    
        [E_test,Y_test,~,~] = SimulateMarkedHawkesMD(end_time, v, zeros(2,1), B, 'const', A);
        E_test1 = cell2mat(E_test(1));
        E_test2 = cell2mat(E_test(2));
        bins = 0:bin_width:end_time;
        bin_sim1 = transpose(histcounts(E_test1, bins));
        bin_sim2 = transpose(histcounts(E_test2, bins)); 
        data = [bin_sim1, bin_sim2];
    end

    disp('Beginning estimation')
    % user set parameters    
    N_monte_carlo = 10;  
    n_times = 10; 
    seed = 0;
    true_data{1} = E_test1; true_data{2} = E_test2;
    n_splits = 10;

    % METHOD 1 (slower but less biased)
    [params_sequential, true_value_est, fval1, grads, hessians, sequential_timer] = MCEM_algorithm_full(data, true_data, N_monte_carlo, n_times, seed, end_time, bin_width, n_splits, 'sequential');
    % METHOD 2 (quicker but more biased)
    [params_uniform, true_value_est, fval1, grads, hessians, uniform_timer] = MCEM_algorithm_full(data, true_data, N_monte_carlo, n_times, seed, end_time, bin_width, n_splits, 'uniform');

    disp('Uniform MCEM Parameter Estimates')
    disp([sum(params_uniform(3:2:end,:),1);sum(params_uniform(4:2:end,:),1)]/n_times)
    disp('Sequential MCEM Parameter Estimates')
    disp([sum(params_sequential(3:2:end,:),1);sum(params_sequential(4:2:end,:),1)]/n_times)


end
