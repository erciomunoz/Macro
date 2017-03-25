%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Figure 4.4 and 4.5 (ABCS of RBCs)
%                            Ercio Munoz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

tic % Start timer

% Parameter values of the model
gamma  = 0.5;   % Weight of consumption leisure in the ut. fn.
alpha  = 0.36;  % Elasticity of output w.r.t. capital
beta   = 0.98;  % Discount factor 
delta  = 0.1;   % Depreciation rate

% Criteria to stop the algorithm
tol    = 0.01;  % Maximum tolerance
maxits = 1000;  % Maximum number of iterations

% Form capital grid
x = [.0001:.05:10]'; % Capital grid
N = length(x);  % Number of grid points for capital
n = [0:0.01:1]; % Labor grid
H = length(n);  % Number of grid points for labor
% Create arrays U(NX,NY) (utility), V(NX) (value fn), Y(NX) (control var.)
U     = zeros(N,N,H);
V_0   = zeros(N,1); 
Y     = zeros(N,1);
labor = zeros(N,1);
% Compute the array U(NX,NY), assigning negative
% values to choices of controls that are not feasible.
k           = repmat(x,1,N);
k_t         = k';
k           = repmat(k,1,1,H);
k_t         = repmat(k_t,1,1,H);
h           = NaN(1,1,H);
h(1,1,:)    = n(:);
h           = repmat(h,N,N,1);
c = (k.^alpha).*(h.^(1-alpha)) + (1-delta)*k - k_t; % Compute consumption
c(find(c<=0)) = NaN; % C must be non-negative
c(find(h==1)) = NaN; % C must be non-negative
U = log(c)+gamma*log(1-h); % Compute utility
U(find(isnan(U))) = -inf; % Set non-feasible decision with -infinity

v_lhs = zeros(N,1);
its = 0;   % Set iterations at zero
dif = 100; % Set the difference higher than the tolerance
while dif > tol && its < maxits 
    for NX = 1:N % Loop over the state variable
        V_1       = repmat(V_0',1,1,H);
        [d_0,i_0] = max(U(NX,:,:)+beta*V_1); % Max lhs V given a state NX
        [d_1,i_1] = max(d_0(1,1,:)); 
        labor(NX) = n(i_1);
        Y(NX)     = k_t(1,i_0(1,1,i_1),i_1);    % Decision k' that maximizes
        v_lhs(NX,1) = d_1;
    end
dif = norm(v_lhs-V_0); % Compute the diference between past and new V
its = its+1; % Counting iterations
V_0 = v_lhs;
end

% Print the number of iterations
if its<maxits
    formatSpec = 'The algorithm stopped after %1.0f iterations. ';
    fprintf(formatSpec,its)
else formatSpec = 'The algorithm reached the max. number of iterations. ';
end

toc % Time elapsed

figure(1)
subplot(1,3,1);
plot(x, V_0,'--*');
title('Converged value function');
xlabel('Capital stock');
ylabel('Value function');

subplot(1,3,2);
plot(x, Y,'--*');
hold on
plot(x, x,'-');
title('Decision rules for capital');
xlabel('Capital stock');
ylabel('Capital stock next period');
legend('Decision rule','45 degree line');

subplot(1,3,3);
plot(x, labor,'--*');
title('Decision rules for labor');
xlabel('Capital stock');
ylabel('Labor decision');
%print figure1.eps