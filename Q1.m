%% Part a

clc; clear;

fprintf("Part A: X^2 calculations\n");

data = load("Q1_data.mat").data;

avg_y = mean(data(:,2:end),2);
std_y = std(data(:,2:end),0,2)/sqrt(length(data(1,2:end)));

x_avg_std_tab = [data(:,1), avg_y, std_y];

disp(x_avg_std_tab)

% power-law relation
pwr_law = @(p,x) p(1).*x.^(p(2));

% chi squared function for power-law fit
chisq_power = @(p) sum( ( (avg_y - pwr_law(p, data(:,1)))./(std_y) ).^2);

% initial parameters (assume direct proportionality)
p_0 = [1, 1];

pmin = fminsearch(chisq_power, p_0);

minX2 = chisq_power(pmin);

N = length(data(:,1));
nc = length(pmin);

% calculate parameter fit errors with x^2+1 method
x2p1 = [];
r_t = [];

for t = 0:0.02:2*pi
    % use r and t to create ellipse of X^2_min + 1
    % define new function in terms of r and t
    chisq_min = @(r) sum(((avg_y - pwr_law([pmin(1)+r*cos(t) pmin(2)+r*sin(t)],data(:,1)))./(std_y)).^2)-(minX2+1);
    r_start = 0;
    r_zero = -1;
    while(r_zero < 0)
        r_start = r_start+2;
        r_zero = fzero(chisq_min, r_start);
    end
    r_t = [r_t; [r_zero*cos(t) r_zero*sin(t)]];
    A = pmin(1)+r_zero*cos(t);
    p = pmin(2)+r_zero*sin(t);
    x2p1 = [x2p1; [A, p]];
end

figure(1)
plot(x2p1(:,1),x2p1(:,2), 'bx')
hold on
plot(pmin(1), pmin(2), 'r+')
xlabel("A")
ylabel("p")
title("X^2_{min} + 1")

A_sigma = (max(x2p1(:,1)) - min(x2p1(:,1)))/2;
p_sigma = (max(x2p1(:,2)) - min(x2p1(:,2)))/2;

[A_mesh, p_mesh] = meshgrid([5:0.1:20],[0:0.01:2]);

for i = 1:length(A_mesh(1,:))
   for j = 1:length(p_mesh(:,1))
      Z(j,i) = chisq_power([A_mesh(j,i) p_mesh(j,i)]); 
   end
end

figure(2)
surf(A_mesh,p_mesh,Z)

% calculate R^2
y_fit = pwr_law(pmin,data(:,1));
R2 = 1 - sum((avg_y - y_fit).^2)/sum((avg_y - ones(length(avg_y),1)*mean(avg_y)).^2);

% print calculation results
fprintf("minimum Chi square of %.3f\n", minX2);
fprintf("minimum Chi square reduced of %.3f\n", minX2/(N-nc));
fprintf("fit parameters A = %.3f +- %.3f, p = %.3f +- %.3f\n",pmin(1), A_sigma, pmin(2), p_sigma);
fprintf("R^2 = %.3f\n",R2);

% calculate best fit points
x_fit = 0.1:0.01:2;
y_fit = pwr_law(pmin,x_fit);

% plot best fit line against data
figure(3)
errorbar(data(:,1),avg_y, std_y, 'k*', 'DisplayName', 'Raw Data')
hold on;
plot(x_fit,y_fit,'b-', 'DisplayName', 'X^2 Fit')

%% Part b

% analytical OLS
fprintf("Analytical OLS Fitting\n");

% set linearized variables
x_tilde = log(data(:,1));
y_tilde = log(avg_y);

% calculate X and Y
X = [(x_tilde.^0) (x_tilde.^1)];
Y = y_tilde;

alpha = X'*X;
alpha_inv = inv(alpha);

param_ols = alpha_inv*X'*Y;

p_ols = param_ols(2);
p_std_ols = sqrt(alpha_inv(2,2));

A_ols = exp(param_ols(1));
A_std_ols = A_ols*sqrt(alpha_inv(1,1));

% calculate R^2
y_fit = pwr_law([A_ols p_ols],data(:,1));
R2 = 1 - sum((avg_y - y_fit).^2)/sum((avg_y - ones(length(avg_y),1)*mean(avg_y)).^2);

% print results
fprintf("The fit parameters were A = %.3f +- %.3f and p = %.3f +- %.3f\n", A_ols, A_std_ols, p_ols, p_std_ols) 
fprintf("R^2 = %.3f\n", R2);

% calculate best fit points
x_fit = 0.1:0.01:2;
y_fit = pwr_law([A_ols p_ols],x_fit);

% plot best fit line against data
plot(x_fit,y_fit,'r-', 'DisplayName', 'OLS Fit')


%% Part c

% analytical WLS
fprintf("Analytical WLS Fitting\n");

% correct uncertainties for log
y_var_tilde = (std_y./avg_y).^2;

% calculate X and Y
X = [(x_tilde.^0)./y_var_tilde (x_tilde.^1)./y_var_tilde];
Y = y_tilde./y_var_tilde;

alpha = X'*X;
alpha_inv = inv(alpha);

param_wls = alpha_inv*X'*Y;

p_wls = param_wls(2);
p_std_wls = sqrt(alpha_inv(2,2));

A_wls = exp(param_wls(1));
A_std_wls = A_wls*sqrt(alpha_inv(1,1));

% calculate R^2
y_fit = pwr_law([A_wls p_wls],data(:,1));
R2 = 1 - sum((avg_y - y_fit).^2)/sum((avg_y - ones(length(avg_y),1)*mean(avg_y)).^2);

% print results
fprintf("The fit parameters were A = %.3f +- %.3f and p = %.3f +- %.3f\n", A_wls, A_std_wls, p_wls, p_std_wls) 
fprintf("R^2 = %.3f\n", R2);

% calculate best fit points
x_fit = 0.1:0.01:2;
y_fit = pwr_law([A_wls p_wls],x_fit);

% plot best fit line against data
plot(x_fit,y_fit,'g-', 'DisplayName', 'WLS Fit')

%% Part d

% label graph
title('Power-Law Curve Fits')
xlabel('x')
ylabel('y')
set(gca, 'XScale', 'log', 'YScale', 'log')
grid on
legend()



