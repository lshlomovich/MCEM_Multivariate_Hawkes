% MCEM Function MODIFIED: NO CLS COMPUTATION

%%% In this version, we allow sampling to be done either sequentially or via uniform sampling

%MC_EM for hawkes process param estimation with parameters:
%v (baseline),alpha (excitation),beta (decay)

function [params, true_value_est, fval1, grads, hessians, mcem_timer] = MCEM_algorithm_full(data, true_data, N_monte_carlo, n_times, seed, end_time, bin_width, n_splits, sampling_method)

    % data:               the count data to be modelled, number of columns
    % should be the dimensionality of the process
    % N_monte_carlo:      number of monte carlo samples
    % n_times:            a maximum number of times to repeat the EM
    % seed:               if wishing to set the seed - 0 if not
    % n_splits:           number of times to try splitting N^{p+1} before
    % selecting the 'best' N^{1}, ..., N^{p}
    

    
    p = size(data,2);
    params = zeros(p*n_times,p); 
    log_likelihood = ones(1,N_monte_carlo); 
    if p==1
        T_true = zeros(N_monte_carlo,sum(data));
    else
        T_true = {};
    end
    weights = log_likelihood; 

    if seed~= 0
        rng(seed)
    end
    
    if p==1
        params(1,:) = sort(rand(1,p));
    elseif p==2
        % note: can make this random
        params = [[0.1;0.1],[0.4,0.2;0.2,0.4],[1.2,0.6;0.6,1.2]];
    elseif p==3
%         params = [[0.1;0.1;0.1],[0.4,0.2,0.1;0.2,0.4,0.1;0.1,0.1,0.4],[1.2,0.6,0.6;0.6,1.2,0.6;0.6,0.6,1.2]];
        init_alpha = rand(3,3); init_beta = init_alpha + rand(3,3);
        while eigs(init_alpha./init_beta, 1, 'lm')-1 >= 0
            init_alpha = rand(3,3); init_beta = init_alpha + rand(3,3);
        end
        params = [rand(3,1), init_alpha, init_beta];
        %params = [[0.1;0.1;0.1]*(0.5*rand(1)+0.5),[0.4,0.2,0.1;0.2,0.4,0.1;0.1,0.1,0.4]*(0.7*rand(1)+0.5),[1.2,0.6,0.6;0.6,1.2,0.6;0.6,0.6,1.2]*(0.7*rand(1)+0.5)];
    end
    
    if iscell(true_data) 
        if p==2 || p==3
            if p==2
                % note: can make this random
                params = [[0.1;0.1],[0.4,0.2;0.2,0.4],[1.2,0.6;0.6,1.2]];
                E_test1 = true_data{1}; E_test2 = true_data{2};
                T_true{1} = E_test1; T_true{2} = E_test2;
            elseif p==3
                params = [[0.1;0.1;0.1],[0.4,0.2,0.1;0.2,0.4,0.1;0.1,0.1,0.4],[1.2,0.6,0.6;0.6,1.2,0.6;0.6,0.6,1.2]];
                E_test1 = true_data{1}; E_test2 = true_data{2}; E_test3 = true_data{3};
                T_true{1} = E_test1; T_true{2} = E_test2; T_true{3} = E_test3;
            end
            
            options = optimoptions('fmincon','MaxIter',2e4,'MaxFunEvals',2e4,'Algo','trust-region-reflective','TolFun',1e-5,'TolX',1e-5,'GradObj','on','Hessian','on','Display','off','DerivativeCheck','off');
            myfunMv_true2 = @(z) complete_likelihood(T_true, z(:,1),z(:,2:2+p-1),z(:,2+p:end), end_time, p);
            tic; [true_value_est,fval,exitflag,output,lambda,grad,hessian] = fmincon(myfunMv_true2, params,[],[],[],[],zeros(p, 2*p+1),[],[],options); toc
            disp('MLE of raw simulated times')
            disp(true_value_est)
            
            options = optimoptions('fmincon','MaxIter',2e4,'MaxFunEvals',2e4,'Algo','interior-point','TolFun',1e-4,'TolX',1e-4,'GradObj','off','Hessian','off','Display','off','DerivativeCheck','off');
            E = sort([E_test1; E_test2]);
            myfunMv_true = @(z) complete_likelihood(E', z(1),z(2),z(3), end_time, 1);
            tic; [super_params,fval,exitflag,output,lambda,grad,hessian] = fmincon(myfunMv_true, [0.6,sum(sum(params(:,2:3))),mean(mean(params(:,4:5)))],[],[],[],[],[],[],[],options); toc
            disp('MLE of raw simulated times')
            disp(super_params)
        end
    else
        true_value_est = 0;
    end

    grads = [];
    hessians = [];

    start_timer = tic;
    
    j = 1;
    eps = 1;
    while j<n_times && eps>0.01
%     for j=1:n_times %set some tolerance and do while loop here
        clear T
        clear log_densities
        disp(j)
        %EXPECTATION STEP 
        weights_raw = zeros(N_monte_carlo,1);
        for i=1:N_monte_carlo
            if mod(i,10)==0
                disp(i)
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Uniform method %
            if strcmp(sampling_method, 'uniform')
                if p==1
                    [T{i}, log_density] = generate_uniform_times(data,p,1);  
                    log_prob_split = 0;  % no splitting necessary 
                else
                    % Simulate an N^(3) and then split to form each of the 2 processes
                    current_params = params(j*p-p+1:j*p,:);

                    baseline_est = current_params(:,1);
                    excitation_est = current_params(:,2:2+p-1);
                    decay_est = current_params(:,2+p:end);
                    summed_data = sum(data,2); % the superposition process
                    [Tsuperpose, log_density(i)] = generate_uniform_times(summed_data,1,bin_width); Tsuperpose = Tsuperpose';  
                end
            elseif strcmp(sampling_method, 'sequential')
                % Simulate an N^(3) and then split to form each of the 2 processes
                current_params = params(j*p-p+1:j*p,:);

                baseline_est = current_params(:,1);
                excitation_est = current_params(:,2:2+p-1);
                decay_est = current_params(:,2+p:end);
                
                % Reparameterise
                baseline_sum = sum(baseline_est);
                decay_mean = mean(mean(decay_est));
                sum_params = [baseline_sum, (1 - (bin_width)*(baseline_sum/sum(mean(data))))*decay_mean, decay_mean];
                
                summed_data = sum(data,2); % the superposition process
                if p > 1
                    [Tsuperpose, ~, log_density_value] = disc_time_hp_grid(summed_data, 1, sum_params, bin_width, 1); 
                    log_density(i) = log_density_value;
                else
                    disp('have not coded this option')
                    return
                end
            end
        
            % Then split T3
            if p==2
                
                [T1, T2, ~] = split_times2_iter(Tsuperpose, data, bin_width, n_splits, current_params, end_time);
                log_prob_split = 0;
                for ii = 1:length(summed_data)
                    val2 = 1/(nchoosek(summed_data(ii),data(ii,1)));
                    log_prob_split = log_prob_split + log(val2);
                end
                T{i}{1} = T1;
                T{i}{2} = T2;
            elseif p==3
                [T1, T2, T3] = split_times3(Tsuperpose, data, bin_width);
                log_prob_split = 0;
                for ii = 1:length(summed_data)
                    val_T1 = 1/(nchoosek(summed_data(ii),data(ii,1)));
                    val_T2 = 1/(nchoosek(summed_data(ii)-data(ii,1),data(ii,2)));
                    log_prob_split = log_prob_split + log(val_T1) + log(val_T2);
                end
                T{i}{1} = T1;
                T{i}{2} = T2;
                T{i}{3} = T3;
            end

            weights_raw(i) = -complete_likelihood(T{i},baseline_est,excitation_est,decay_est,end_time,p)-(log_density(i));   % +/-end time?
           
        end
        %% Normalise weights: mean(weights) is the IS estimator for the messy integral (normalisation constant)

        if strcmp(sampling_method, 'uniform')
            % For this method weights need to be rescaled to not be Inf
            weights_og = weights_raw;
            weights_raw = weights_og-min(weights_og); 
            weights_raw = exp(weights_raw);
            weights_raw_rescaled = weights_raw./sum(weights_raw);
        elseif strcmp(sampling_method, 'sequential')
            weights_og = weights_raw;
            weights_raw = weights_og-min(weights_og); 
            weights_raw = exp(weights_raw);
            weights_raw_rescaled = weights_raw./sum(weights_raw);
        end
        

        %MAXIMISATION STEP
        if p==2 || p==3
%             zhatMv = params(end-p+1:end,:);
            if size(params,1) == 2
                zhatMv = params(end-p+1:end,:);
            elseif p == 2
                zhatMv = [mean(params(3:2:end,:),1); mean(params(4:2:end,:),1)];
            else
                zhatMv = params(end-p+1:end,:);
            end

            options1 = optimoptions('fmincon','MaxIter',1e3,'MaxFunEvals',1e3,'Algo','trust-region-reflective','TolFun',1e-5,'TolX',1e-4,'GradObj','on','Hessian','on','Display','off','DerivativeCheck','off');

            myfunMv_alt = @(z) complete_sum_likelihood(T, z(:,1),z(:,2:2+p-1),z(:,2+p:end), end_time, weights_raw_rescaled,p, N_monte_carlo); % weights3 doesnt have NaNs, weights_alt most logical, weights_rescaled is an attempt to stop NaNs
            weights1_time=tic; [new_params, fval1(j)] = fmincon(myfunMv_alt,zhatMv,[],[],[],[],[],[],[],options1); w1_time = toc(weights1_time); disp(w1_time)
            % Test run time without using the hessian
            %weights1_time=tic; [new_params_alt, fval1_alt(j)] = fmincon(myfunMv_alt,zhatMv,[],[],[],[],[],[],[],options1_nohess); w1_time = toc(weights1_time); disp(w1_time)
            disp(new_params)
            eps = norm(new_params - params(end-1:end,:));
            params = [params;new_params];
            disp(fval1(j))
            %[sum(params(3:2:end,:),1);sum(params(4:2:end,:),1)]/j
        else
            zhat = params(j,:); 
            myfun = @(z) complete_sum_likelihood(T, z(:,1),z(:,2:2+p-1),z(:,2+p:end), end_time, weights_raw_rescaled,p, N_monte_carlo); 

            options1 = optimoptions('fmincon','MaxIter',2e4,'MaxFunEvals',2e4,'Algo','trust-region-reflective','TolFun',1e-4,'TolX',1e-4,'GradObj','on','Hessian','on','Display','off','DerivativeCheck','off');

            weights1_time=tic; [params(j+1,:), fval1(j),~,~,~,grad,hessian] = fmincon(myfun,zhat,[],[],[],[],[],[],[],options1); w1_time = toc(weights1_time); disp(w1_time)
            grads = [grads;grad];
            hessians = [hessians;hessian];
            disp(params(j+1,:))
        end
        
        disp('eps is')
        disp(eps)
        j = j+1;
    end
     
    mcem_timer = toc(start_timer);
    
end

%% Functions - MV simulation ideas
function [T1, T2] = split_times2(T3, data, bin_width)
% For the bivariate case 
T1 = [];
for i=1:length(data)
        if data(i,1) ~= 0
                range_times_ub = i*bin_width; 
                range_times_lb = range_times_ub - bin_width;
                sample_times = T3(T3 < range_times_ub);
                sample_times = sample_times(sample_times >= range_times_lb);

                T1 = [T1;sort(datasample(sample_times,data(i,1),'Replace',false))];
%                 diffs = diff(sample_times);
%                 [~, max_idxs] = maxk(diffs,data(i,1));
%                 T1 = [T1;sort(sample_times(max_idxs+1))];
        end
end
T2 = setdiff(T3,T1);
end

function [T1, T2, log_lls] = split_times2_iter(T3, data, bin_width, N, params, end_time)
% Iterated splitting and selecting of the best one
% For the bivariate case 
counter=1;
log_lls = zeros(N,1); % currently unused
p = size(params,1);
while counter<=N

    T1 = [];
    for i=1:length(data)
            if data(i,1) ~= 0
                    range_times_ub = i*bin_width; 
                    range_times_lb = range_times_ub - bin_width;
                    sample_times = T3(T3 < range_times_ub);
                    sample_times = sample_times(sample_times >= range_times_lb);

                    T1 = [T1;sort(datasample(sample_times,data(i,1),'Replace',false))];
    %                 diffs = diff(sample_times);
    %                 [~, max_idxs] = maxk(diffs,data(i,1));
    %                 T1 = [T1;sort(sample_times(max_idxs+1))];
            end
    end
    T2 = setdiff(T3,T1);

    T{counter}{1} = T1;
    T{counter}{2} = T2;
    counter = counter+1;
end

for counter_idx = 1:N
    log_lls(counter_idx) = complete_likelihood(T{counter_idx}, params(:,1),params(:,2:2+p-1),params(:,2+p:end), end_time,p);
end

[~,idx] = min(log_lls);

T1 = T{idx}{1};
T2 = T{idx}{2};

end

function [T1, T2, T3] = split_times3(T4, data, bin_width)
% For the trivariate case (will eventually make general)
T1 = [];
for i=1:length(data)
        % Assign T1 first
        if data(i,1) ~= 0
                range_times_ub = i*bin_width; 
                range_times_lb = range_times_ub - bin_width;
                sample_times = T4(T4 < range_times_ub);
                sample_times = sample_times(sample_times >= range_times_lb);
                % the times to be added to T1 are:
                times_added = sort(datasample(sample_times,data(i,1),'Replace',false));
                T1 = [T1;times_added];
        end
end
T5 = setdiff(T4,T1);
% Now sample T2
T2 = [];
for i=1:length(data)
        % Assign T2
        if data(i,2) ~= 0
                range_times_ub = i*bin_width; 
                range_times_lb = range_times_ub - bin_width;
                sample_times = T5(T5 < range_times_ub);
                sample_times = sample_times(sample_times >= range_times_lb);
                % the times to be added to T1 are:
                times_added = sort(datasample(sample_times,data(i,2),'Replace',false));
                T2 = [T2;times_added];
        end
end
T3 = setdiff(T5,T2);
end
