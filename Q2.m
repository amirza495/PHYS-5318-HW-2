clc; clear;

fprintf("Part A: X^2 calculations\n");

data = load("Q2_data.mat").data;

avg_y = mean(data(:,2:end),2);
std_y = std(data(:,2:end),0,2)/sqrt(length(data(1,2:end)));

x_avg_std_tab = [data(:,1), avg_y, std_y];

disp(x_avg_std_tab)

% dna changes relation
dna_change = @(p,x) (p(1) - (p(1) - 0.355).*exp(-p(2).*x));

% chi squared function for power-law fit
chisq_power = @(p) sum(((avg_y - dna_change(p,data(:,1)))./(std_y)).^2);

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
    chisq_min = @(r) sum(((avg_y - dna_change([pmin(1)+r*cos(t) pmin(2)+r*sin(t)],data(:,1)))./(std_y)).^2)-(minX2+1);
    r_start = 0;
    r_zero = -0.00001;
    while(r_zero < 0)
        r_start = r_start+0.001;
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

[A_mesh, p_mesh] = meshgrid([0.399:0.0001:0.402],[0.05:0.0001:0.4]);

for i = 1:length(A_mesh(1,:))
   for j = 1:length(p_mesh(:,1))
      Z(j,i) = chisq_power([A_mesh(j,i) p_mesh(j,i)]); 
   end
end

figure(2)
surf(A_mesh,p_mesh,Z)

% calculate R^2
y_fit = dna_change(pmin,data(:,1));
R2 = 1 - sum((avg_y - y_fit).^2)/sum((avg_y - ones(length(avg_y),1)*mean(avg_y)).^2);

% print calculation results
fprintf("minimum Chi square of %.3f\n", minX2);
fprintf("minimum Chi square reduced of %.3f\n", minX2/(N-nc));
fprintf("fit parameters b_eq = %.4f +- %.6f, k = %.3f +- %.3f\n",pmin(1), A_sigma, pmin(2), p_sigma);
fprintf("R^2 = %.3f\n",R2);

% calculate best fit points
x_fit = 0.1:0.01:max(data(:,1));
y_fit = dna_change(pmin,x_fit);

% plot best fit line against data
figure(3)
errorbar(data(:,1),avg_y, std_y, 'k*', 'DisplayName', 'Raw Data')
hold on;
plot(x_fit,y_fit,'b-', 'DisplayName', 'X^2 Fit')
xlabel("time (min)")
ylabel("extension per base pair")
title("Extension per base pair vs time")