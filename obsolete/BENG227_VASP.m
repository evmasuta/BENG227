%% BENG 227 MIDTERM PROJECT Evan Masutani
%% HOUSEKEEPING
clear all
close all
clc
%% SIMULATION RUN PARAMETERS
% step sizes
dr = 1;
% Keep isotropic to allow for downstream 2D Laplacian filter
dc = dr;
% Declare an isotropic spacer, dx
dx = dr;
dt = 0.01;
% limits
R_dim = 40;
C_dim = 10;
T_dim = 100;
% N, +2 denotes phantom points to help with Neumann B.C.'s
N_r = round(R_dim / dr) + 2;
N_c = round(C_dim / dc) + 2;
N_t = round(T_dim / dt);
%% CONSTANTS
R = 0.14;
delta = 5;
E = 0.2;
eta_M = 1;
eta_B = 1;
K = 1;
k_2 = 1;
k_3 = 10;
B0_c = 3.6;
D_A = 3;
D= 0.2;
epsilon = 0.1;

%% INITIAL CONDITIONS
A_b = 2 * ones(N_r,N_t);
A_c = 5 * ones(N_r,N_c,N_t);
% Constant source; set both phantom and real boundary to delta
A_c(:,N_c,:) = delta;
A_c(:,N_c - 1,:) = delta;
A_m = 1*ones(N_r,N_c,N_t);
V = 0*ones(N_r,N_t);
B = 1.5*B0_c*ones(N_r,N_t);
M = 1*ones(N_r,N_c,N_t);
noise = 0.01;
%% DEFINE ddc and ddr matrices; implicitly have Neumann programmed in
ddc = zeros(N_c,N_c);
for i=2:1:N_c-1
   ddc(i, i-1) = -1;
   ddc(i, i+1) = 1;
end

ddr = zeros(N_r,N_r);
for i=2:1:N_r-1
    ddr(i,i-1) = -1;
    ddr(i,i+1) = 1;
end
%% RUN WITH BOUNDARY CONDITIONS
for t=1:1:N_t - 1
    t
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
%     G = ones(N_r,1) + (A_b(:,t) + A_c(:,2,t) + A_m(:,2,t)) .* B(:,t) ./ (ones(N_r,1) + M(:,2,t) + K * B(:,t));
    G = ones(N_r,1) + (A_b(:,t)) .* B(:,t) ./ (ones(N_r,1) + M(:,2,t) + K * B(:,t));
    %% CALCULATE B
    % No flux BC
    B(1,t) = B(3,t);
    B(N_r,t) = B(N_r - 2,t);
    % No need to transpose since natively row vector
    dBdx = ddr * B(:,t) / (2*dx);
    gdBdx = 1./G .* dBdx;
    diffusion_B = epsilon^2 * ddr * gdBdx / (2*dx);
    accumulation_B = ones(N_r,1) + eta_B * V(:,t);
    loss_B = -B(:,t) ./ (ones(N_r,1) + A_b(:,t));
    noise_B = 0.1 * 2 * (rand(N_r,1) - ones(N_r,1));
    B(:,t+1) = B(:,t) + dt/epsilon * (diffusion_B + accumulation_B + loss_B);

    %% CALCULATE A_b
    A_b(:,t+1) = A_b(:,t) + dt * k_2 * (B(:,t) .* A_c(:,2,t) - A_b(:,t));
    %% CALCULATE A_c
    % Calculate phantom points for A_c; Neumann B.C. due to exchange with
    % leading edge
    A_c(:,1,t) = A_c(:,3,t) - 2*dx/D_A * k_2 * (B(:,t) .* A_c(:,2,t) - A_b(:,t));
%     A_c(:,1,t) = A_c(:,3,t);
    % Dirchlet B.C.
    A_c(:,N_c,t) = delta;
    A_c(:,N_c-1,t) = delta;
    % Physical edge Neumann B.C.
    A_c(1,:,t) = A_c(3,:,t);
    A_c(N_r,:,t) = A_c(N_r - 2,:,t);

    % Calculate second derivative Laplacian diffusion term
    diffusion_A_c = D_A * del2(A_c(:,:,t),dx);
    % Calculate reaction terms
    reaction_A_c = k_3 * (A_m(:,:,t) - K * M(:,:,t) .* A_c(:,:,t));
    % Calculate A_c at t+dt
    A_c(:,:,t+1) = A_c(:,:,t) + dt * (diffusion_A_c + reaction_A_c);

    %% CALCULATE M this part is fucked
    % Set leading edge B.C.; Dirchlet B.C. set phantom = real boundary
    if sum(V(:,t) > 0)
        M(:,2,t) = R * B(:,t) ./ V(:,t);
    end
    M(:,1,t) = M(:,2,t);
    % Set implied B.C.; Neumann at col=N_c
    M(:,N_c,t) = M(:,N_c - 2,t);
    for row = 1:1:N_r
        % lots of tranposes here because it's natively dealing with row vectors
        % while finite difference outputs a column vector
        % signs might make a difference here...
        M(row,:,t+1) = M(row,:,t) + transpose(dt * (V(row,t) * 1 / (2* dx) * ddc * transpose(M(row,:,t))));
    end
    %% CALCULATE A_m
    % Set leading edge B.C.; Dirchlet B.C. set phantom = real boundary
    A_m(:,1,t) = 0;
    A_m(:,2,t) = 0;
    % Set implied B.C.; Neumann at col=N_c
    A_m(:,N_c,t) = A_m(:,N_c - 2,t);
    for row = 1:1:N_r
        % lots of tranposes here because it's natively dealing with row vectors
        % while finite difference outputs a column vector
        A_m(row,:,t+1) = A_m(row,:,t) + dt*(-k_3 * (A_m(row,:,t) - K .* M(row,:,t) .* A_c(row,:,t)) + transpose(V(row,t) * 1 / (2* dx) * ddc * transpose(A_m(row,:,t))));
    end
end

figure
tvec = 1:1:N_t-1;
plot(tvec,V(10,1:N_t-1))
title('V')
figure
plot(tvec,B(10,1:N_t-1))
title('B')
figure
plot(tvec,A_b(10,1:N_t-1))
title('A_b')
figure
plot(tvec,squeeze(A_c(10,3,1:N_t-1)))
title('A_c')
figure
plot(tvec,squeeze(A_m(10,3,1:N_t-1)))
title('A_m')
figure
plot(tvec,squeeze(M(10,3,1:N_t-1)))
title('M')
% figure
% plot([2:1:N_t-1],V(10,2:2))
