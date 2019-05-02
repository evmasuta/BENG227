%% BENG 227 MIDTERM PROJECT Evan Masutani
%% HOUSEKEEPING
clear all
close all
clc
%% SIMULATION RUN PARAMETERS
% step sizes
dr = 0.1;
% Declare an isotropic spacer, dx
dx = dr;
dt = 0.01;
% limits
R_dim = 40;
T_dim = 500;
% N, +2 denotes phantom points to help with Neumann B.C.'s
N_r = round(R_dim / dr) + 2;
N_t = round(T_dim / dt);
%% CONSTANTS
R = 0.28;
delta = 0.9;
epsilon = 0.1;
K = 1;
eta_B = 1;
% eta_M and eta_A not clearly labeled
eta_M = 1;
eta_A = 0.2;
E = 0.1;
B0_c = 4;
theta = 0.05;
noise = 0.01;
%% INITIAL CONDITIONS
B = 0 * ones(N_r,N_t);
A = 0 * ones(N_r,N_t);
M = 0 * ones(N_r,N_t);
V = 0 * ones(N_r,N_t);

% Gradient I.C.
for i = 1:1:N_r
   B(i,1) = 20 - i/N_r * 20;
   A(i,1) = 20 - i/N_r * 20;
end
% B(1:100,1) = 10;
% A(1:100,1) = 10;
% M(1:100,1) = 0;
%% DEFINE ddr matrices; implicitly have Neumann programmed in

ddr = zeros(N_r,N_r);
for i=2:1:N_r-1
    ddr(i,i-1) = -1;
    ddr(i,i+1) = 1;
end
ddr = 1 / (2*dx) * ddr;


%% RUN WITH BOUNDARY CONDITIONS
for t=1:1:N_t - 1
    %% CALCULATE V
    if sum(V(:,t)) > 0
        B_c = B0_c * (1 + E);
    else
        B_c = B0_c;
    end
    V(:,t) = ones(N_r,1) - (B_c * ones(N_r,1) ./ B(:,t)).^8;
    for vstep = 1:1:N_r
       if V(vstep,t) < 0
          V(vstep,t) = 0; 
       end
    end
    % G is kinda weird, consider changing if fails
    G = ones(N_r,1) + (A(:,t) .* B(:,t)) ./ (ones(N_r,1) + M(:,t) + K * B(:,t));
    %% CALCULATE B
    % No flux BC
    B(1,t) = B(3,t);
    B(N_r,t) = B(N_r - 2,t);
    % No need to transpose since natively row vector
    dBdx = ddr * B(:,t);
    gdBdx = 1./G .* dBdx;
    diffusion_B = epsilon^2 * ddr * gdBdx;
    accumulation_B = ones(N_r,1) + eta_B * V(:,t);
    loss_B = -B(:,t) ./ (G);
    noise_B = 0.1 * 2 * (rand(N_r,1) - ones(N_r,1));
    B(:,t+1) = B(:,t) + dt/epsilon * (diffusion_B + accumulation_B + loss_B + noise_B);

    %% CALCULATE A
    accumulation_A = ones(N_r,1) * delta;
    loss_A = -1 ./ (ones(N_r,1) + M(:,t) + K .* B(:,t)) .* (ones(N_r,1) + eta_A .* V(:,t) + eta_M .* M(:,t) .* V(:,t)) .* A(:,t);
    A(:,t+1) = A(:,t) + dt * (accumulation_A + loss_A);
    %% CALCULATE M
    accumulation_M = R * B(:,t);
    loss_M = -1 * (ones(N_r,1) * theta + eta_M * V(:,t)) .* M(:,t);
    M(:,t+1) = M(:,t) + dt * (accumulation_M + loss_M);
end

tvec = 1:1:N_t;
% figure
% plot(tvec,M(10,:))
% figure
% plot(tvec,A(10,:))
% figure
% plot(tvec,B(10,:))
figure
imagesc(V(:,1:N_t-1))
colorbar
figure
imagesc(B(:,1:N_t-1))
colorbar
% figure
% plot(tvec,V(10,:));
figure
imagesc(A(:,1:N_t-1))
colorbar
figure
imagesc(M(:,1:N_t-1))
colorbar
figure
plot(tvec,B(10,:))
figure
V_scale = V ./ max(max(V));
imshow(V)
