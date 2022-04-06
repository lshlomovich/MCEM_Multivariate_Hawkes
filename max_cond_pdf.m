function [pdf_val,interim_pdf_val] = max_cond_pdf(miss_times, history, bin_width, max_bin, v, alpha, beta)
% This function computes the likelihood of the current bin given the
% history of the process. That is, with the parameters v,alpha,beta, and
% the history of the process, it can take proposed missing times and
% compute their likelihood. Thus, maximising this function allows us to
% resdistribute points. Note that we truncate this likelihood to ensure it
% sits between the bin bounds.

% m events in bin to be simulated
% miss_times is row vector of missing time points 
% needs adjusting for any bin width TO DO

% if ~isempty(future)
%     final_bin = ceil(future(end));
% else
%     final_bin = max_bin;
% end

m = length(miss_times);
if any(sort(miss_times)~=miss_times)
    pdf_val = inf;
    return
end

% if min(diff(miss_times)) < 0.01
%     pdf_val = inf;
%     return
% end

if (min(miss_times)<=(max_bin-bin_width) || max(miss_times)>=max_bin)
    pdf_val = inf;
    return
end

if length(unique(miss_times))~=m
    pdf_val = inf;
    return
end

all_events = [0;history;miss_times'];
%all_events = [0;history;miss_times'; future]; % if the zero is removed, adjust point_n

point_n = length(history)+1;
product_term = 1;
for i=1:m %+length(future)   %m   % length(all_events)-length(history)-1
    first_product = v + alpha*sum(exp(-beta*(all_events(point_n+i)-all_events(2:point_n+i-1))));
    second_product = exp(alpha/beta * sum( exp(-beta*(all_events(point_n+i)-all_events(2:point_n+i-1))) - exp(-beta*(all_events(point_n+i-1)-all_events(2:point_n+i-1)))));
    product_term = product_term * first_product * second_product;
end
pdf_val = product_term*exp(-v*(all_events(end)-all_events(point_n)));
interim_pdf_val = pdf_val;
%test_prod_ub = 1 - exp(-v*(next_bin-all_events(end)))*exp(alpha/beta * (sum(exp(-beta*(next_bin-all_events(2:end)))) - sum(exp(-beta*(all_events(end)-all_events(2:end))))));
%pdf_val_test = pdf_val*(1-test_prod_ub);
%interim_pdf_val = -pdf_val*(1-test_prod_ub);

% TRUNCATION METHOD - uncomment
upper_cdf = 1;
%ub_times = [0;history;transpose(miss_times(2:end));max_bin]; % OG
%ub_times = [0;history;transpose(miss_times(2:end));max_bin;future]; %OG % +1 to deal with  edge effect?
ub_times = [0;history;transpose(miss_times(1:end-1));max_bin];
middle_prod = 1;
for i=1:m
    %ith_prod = 1 - exp(-v*(ub_times(end-i+1)-ub_times(end-i)))*exp(alpha/beta * (sum(exp(-beta*(ub_times(end-i+1)-ub_times(1:end-i)))) - sum(exp(-beta*(ub_times(end-i)-ub_times(1:end-i)))))); %OG
    ith_prod = 1 - exp(-v*(ub_times(point_n+i)-ub_times(point_n+i-1)))*exp(alpha/beta * (sum(exp(-beta*(ub_times(point_n+i)-ub_times(2:point_n+i-1)))) - sum(exp(-beta*(ub_times(point_n+i-1)-ub_times(2:point_n+i-1)))))); % same as OG
    % OR
    %ith_prod = 1 - exp(-v*(ub_times(point_n+i)-all_events(point_n+i-1)))*exp(alpha/beta * (sum(exp(-beta*(ub_times(point_n+i)-all_events(2:point_n+i-1)))) - sum(exp(-beta*(all_events(point_n+i-1)-all_events(2:point_n+i-1))))));
    %ith_prod = 1 - exp(-v*(max_bin-all_events(end-i)))*exp(alpha/beta * (sum(exp(-beta*(max_bin-all_events(1:end-i)))) - sum(exp(-beta*(all_events(end-i)-all_events(1:end-i))))));
%     if i==1 % for rectangular idea
%         t1_val = ith_prod;
%     end
%     if i==m % for rectangular idea
%         bplus_val = ith_prod;
%     end
    if (i>=2 && i <= m-1)
        middle_prod = middle_prod * ith_prod;
    end
    upper_cdf = ith_prod * upper_cdf;
end
lb_times = [0;history;transpose(miss_times(1:end))];
final_cdf_val = 1 - exp(-v*(lb_times(point_n+m)-lb_times(point_n+m-1)))*exp(alpha/beta * (sum(exp(-beta*(lb_times(point_n+m)-lb_times(2:point_n+m-1)))) -sum(exp(-beta*(lb_times(point_n+m-1)-lb_times(2:point_n+m-1))))));
first_cdf_val = 1 - exp(-v*(max_bin-bin_width-lb_times(point_n)))*exp(alpha/beta * (sum(exp(-beta*(max_bin-bin_width-lb_times(2:point_n)))) -sum(exp(-beta*(lb_times(point_n)-lb_times(2:point_n))))));
lower_cdf = middle_prod * final_cdf_val * first_cdf_val;


% PREVIOUS TRUNCATION: sometimes neg. - unused now
% lower_cdf = 1;
% lb_times = [0;history;max_bin-1;transpose(miss_times(1:end-1))]; % OG
% %lb_times = [0;history;max_bin-1;transpose(miss_times(1:end-1));future]; %OG
% %lb_times = [0;history;max_bin-1;transpose(miss_times(2:end))];
% 
% for i=1:m  % if this is m-1 then this is much larger than upper_cdf
%     ith_prod = 1 - exp(-v*(lb_times(point_n+i)-lb_times(point_n+i-1)))*exp(alpha/beta * (sum(exp(-beta*(lb_times(point_n+i)-lb_times(2:point_n+i-1)))) -sum(exp(-beta*(lb_times(point_n+i-1)-lb_times(2:point_n+i-1)))))); % OG
%     %ith_prod = 1 - exp(-v*(lb_times(point_n+i)-all_events(point_n+i-1)))*exp(alpha/beta * (sum(exp(-beta*(lb_times(point_n+i)-all_events(2:point_n+i-1)))) -sum(exp(-beta*(all_events(point_n+i-1)-all_events(2:point_n+i-1)))))); 
%     %ith_prod = 1 - exp(-v*(lb_times(end-i+1)-all_events(end-i)))*exp(alpha/beta * (sum(exp(-beta*(lb_times(end-i+1)-all_events(1:end-i)))) - sum(exp(-beta*(all_events(end-i)-all_events(1:end-i))))));
%     %ith_prod = 1 - exp(-v*(max_bin-1-all_events(end-i)))*exp(alpha/beta * (sum(exp(-beta*(max_bin-1-all_events(1:end-i)))) - sum(exp(-beta*(all_events(end-i)-all_events(1:end-i))))));
%     lower_cdf = ith_prod * lower_cdf;
% end

% Testing rectanglar idea:  seems to underestimate parameters, unused now
% lower_cdf is F(b-,...,t_m)
% upper_cdf is F(t_1,...,b+)
% cross1_cdf is F(t_1,...,t_m) = t1_val*middle_prod*final_cdf_val
% cross2_cdf is F(b-,...,b+) = first_cdf_val*middle_prod*bplus_val
% cross1_cdf = t1_val*middle_prod*final_cdf_val;
% cross2_cdf = first_cdf_val*middle_prod*bplus_val;
% total_cdf = upper_cdf + lower_cdf - cross1_cdf - cross2_cdf;

if upper_cdf-lower_cdf ~= 0
    conditional_pdf_val = pdf_val/(upper_cdf-lower_cdf); %OG with abs()   %%%%((upper_cdf-lower_cdf+ref_cdf));  %
else
    disp('0 for max_cond_pdf denom')
    conditional_pdf_val = pdf_val;
end
%conditional_pdf_val = pdf_val*middle_prod/(upper_cdf-lower_cdf); 
%conditional_pdf_val = pdf_val_test/cdf_total_diff;
%conditional_pdf_val = pdf_val/total_cdf; %total_cdf;

% Negate so that getting minimum gives desired result
%pdf_val = -pdf_val;
pdf_val = -conditional_pdf_val;


end




% % Truncate %%%%%%
% upper_cdf = 1;
% upper_cdf_test = 1;
% ub_times = [0;history;transpose(miss_times(2:end));max_bin]; %OG
% %ub_times = [0;history;transpose(miss_times(1:end-1));max_bin];
% A_ub = zeros(length(ub_times),1);
% A_ub(1:point_n) = A_s(1:point_n);
% for tp = point_n:length(ub_times)-1
%     %tp_diff_ub = ub_times(tp+1) - all_events(tp); 
%     tp_diff_ub = ub_times(tp+1) - ub_times(tp); %OG
%     A_ub(tp+1) = (exp(-beta*tp_diff_ub))*(1+A_ub(tp)); %OG
%     %A_ub(tp+1) = (exp(-beta*tp_diff_ub))*(1+A_s(tp)); 
% end
% for i=1:m
%     %prod_term = 1 - exp(-v*(ub_times(point_n+i)-all_events(point_n+i-1)))*exp(alpha/beta * (A_ub(point_n + i)-A_ub(point_n+i-1)-1));
%     prod_term = 1 - exp(-v*(ub_times(point_n+i)-ub_times(point_n+i-1)))*exp(alpha/beta * (A_ub(point_n + i)-A_ub(point_n+i-1)-1)); %OG
%     %test_prod = 1 - exp(-v*(ub_times(end-i+1)-all_events(end-i)))*exp(alpha/beta * (sum(exp(-beta*(ub_times(end-i+1)-all_events(1:end-i)))) - sum(exp(-beta*(all_events(end-i)-all_events(1:end-i))))));
%     upper_cdf = prod_term * upper_cdf;
%     %upper_cdf_test = test_prod * upper_cdf_test;
% end
% % 
% lower_cdf = 1;
% lower_cdf_test = 1;
% lb_times = [0;history;max_bin-1;transpose(miss_times(1:end-1))]; %OG
% %lb_times = [0;history;max_bin-1;transpose(miss_times(2:end))];
% %lb_times = [0;history;transpose(miss_times(1:end))];
% A_lb = zeros(length(lb_times),1);
% A_lb(1:point_n) = A_s(1:point_n);
% for tp = point_n:length(lb_times)-1
%     %tp_diff_lb = lb_times(tp+1) - all_events(tp);
%     tp_diff_lb = lb_times(tp+1) - lb_times(tp); %OG
%     A_lb(tp+1) = (exp(-beta*tp_diff_lb))*(1+A_lb(tp)); %OG
%     %A_lb(tp+1) = (exp(-beta*tp_diff_lb))*(1+A_s(tp));
% end
% for i=1:m
%     %prod_term = 1 - exp(-v*(lb_times(point_n+i)-all_events(point_n+i-1)))*exp(alpha/beta * (A_lb(point_n + i)-A_lb(point_n+i-1)-1));
%     prod_term = 1 - exp(-v*(lb_times(point_n+i)-lb_times(point_n+i-1)))*exp(alpha/beta * (A_lb(point_n + i)-A_lb(point_n+i-1)-1)); %OG
%     %test_prod = 1 - exp(-v*(lb_times(end-i+1)-all_events(end-i)))*exp(alpha/beta * (sum(exp(-beta*(lb_times(end-i+1)-all_events(1:end-i)))) - sum(exp(-beta*(all_events(end-i)-all_events(1:end-i))))));
%     lower_cdf = prod_term * lower_cdf;
%     %lower_cdf_test = test_prod * lower_cdf_test;
% end




% % % Think this cross idea could be used? 
% perm_array = permn(1:2,m);
% % % Remove rows which are all 1s or all 2s
% perm_array = perm_array(2:end-1,:);
% all_cdfs = zeros(length(perm_array),1);
% for i=1:length(perm_array)
%     cdf_vals = zeros(1,m);
%     for j=1:m
%         if perm_array(i,j)==1 % FIX THIS
%             cdf_vals(j) = lb_times(point_n+j); %max(lb_times(point_n+j),cdf_vals(max(j-1,1)));
%         elseif perm_array(i,j)==2
%             cdf_vals(j) = ub_times(point_n+j); %max(ub_times(point_n+j),cdf_vals(max(j-1,1)));
%         end
%     end
%     cross_times = [0;history;cdf_vals']; %% note need to check terms where next t is larger than previous
%     cross_cdf = 1;
%     for jj=1:m
%         ith_prod = 1 - exp(-v*(cross_times(end-jj+1)-cross_times(end-jj)))*exp(alpha/beta * (sum(exp(-beta*(cross_times(end-jj+1)-cross_times(2:end-jj)))) -sum(exp(-beta*(cross_times(end-jj)-cross_times(2:end-jj)))))); % OG
%         cross_cdf = ith_prod * cross_cdf;
%     end
%     all_cdfs(i) = cross_cdf;
% end
% final = upper_cdf+lower_cdf-sum(all_cdfs);
% 
% % method just for 2 events
% cross1 = [0;history;min(ub_times(end-1),lb_times(end));lb_times(end)];
% cross1_cdf = 1;
% for i=1:m
%     ith_prod = 1 - exp(-v*(cross1(end-i+1)-cross1(end-i)))*exp(alpha/beta * (sum(exp(-beta*(cross1(end-i+1)-cross1(2:end-i)))) -sum(exp(-beta*(cross1(end-i)-cross1(2:end-i)))))); % OG
%     cross1_cdf = ith_prod * cross1_cdf;
% end
% cross2 = [0;history;lb_times(end-1);ub_times(end)];
% cross2_cdf = 1;
% for i=1:m
%     ith_prod = 1 - exp(-v*(cross2(end-i+1)-cross2(end-i)))*exp(alpha/beta * (sum(exp(-beta*(cross2(end-i+1)-cross2(2:end-i)))) -sum(exp(-beta*(cross2(end-i)-cross2(2:end-i)))))); % OG
%     cross2_cdf = ith_prod * cross2_cdf;
% end
%conditional_pdf_val = pdf_val/(upper_cdf+lower_cdf-cross1_cdf-cross2_cdf);

% b_times_ref = [0;history;transpose(miss_times(1:end))];
% ref_cdf = 1;
% for i=1:m
%     %ith_prod = 1 - exp(-v*(ub_times(end-i+1)-ub_times(end-i)))*exp(alpha/beta * (sum(exp(-beta*(ub_times(end-i+1)-ub_times(1:end-i)))) - sum(exp(-beta*(ub_times(end-i)-ub_times(1:end-i)))))); %OG
%     ith_prod = 1 - exp(-v*(b_times_ref(point_n+i)-b_times_ref(point_n+i-1)))*exp(alpha/beta * (sum(exp(-beta*(b_times_ref(point_n+i)-b_times_ref(2:point_n+i-1)))) - sum(exp(-beta*(b_times_ref(point_n+i-1)-b_times_ref(2:point_n+i-1)))))); % same as OG
%     % OR
%     %ith_prod = 1 - exp(-v*(ub_times(point_n+i)-all_events(point_n+i-1)))*exp(alpha/beta * (sum(exp(-beta*(ub_times(point_n+i)-all_events(2:point_n+i-1)))) - sum(exp(-beta*(all_events(point_n+i-1)-all_events(2:point_n+i-1))))));
%     %ith_prod = 1 - exp(-v*(max_bin-all_events(end-i)))*exp(alpha/beta * (sum(exp(-beta*(max_bin-all_events(1:end-i)))) - sum(exp(-beta*(all_events(end-i)-all_events(1:end-i))))));
%     ref_cdf = ith_prod * ref_cdf;
% end



