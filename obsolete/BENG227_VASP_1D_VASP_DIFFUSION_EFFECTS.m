%% BENG 227 MIDTERM PROJECT Evan Masutani
%% HOUSEKEEPING
clear all;
close all;
clc;
%% SIMULATION RUN PARAMETERS
% step sizes
dr = 1;
% Declare an isotropic spacer, dx
dx = dr;
dt = 0.1;
% limits; small size for replicating figures
R_dim = 40; %microns
T_dim = 500;
% N, +2 denotes phantom points to help with Neumann B.C.'s
N_r = round(R_dim / dr) + 2;
N_t = round(T_dim / dt);
tvec = 1:1:N_t;
%% CONSTANTS
R = 0.14;
delta = 0.8;
epsilon = 0.1;
K = 1;
eta_B = 1;
% eta_M and eta_A not clearly labeled
eta_M = 1;
eta_A = 1;
E = 0.1;
B0_c = 4;
theta = 0.2;
noise = 0.01;
% Additional diffusion term for VASP; define based on epsilon2 =
% D_A/cell_l^2/gamma
% VASP estimated diffusivity
D_A = 0.02; %um/sec
cell_l = 40; %um
gamma = 1; % sec
% Non-dimensionalized lateral diffusivity
epsilon2 = D_A / cell_l^2 / gamma * 2;
%% INITIAL CONDITIONS; refer to icy.m for usage
BMaxODE = 8.80342361002194;
BMinODE = 2.46024797933793;
AMaxODE = 5.79255834170108;
AMinODE = 2.26789894792708;
MMaxODE = 2.69213920749395;
MMinODE = 0.677418361848779;
B = noise * rand(N_r,N_t);
A = zeros(N_r,N_t);
M = zeros(N_r,N_t);
V = zeros(N_r,N_t);

% INITIAL CONDITIONS SENT BY AUTHORS
% DON'T YOU DARE TOUCH THESE I.C.s FOR BAM OR I'LL BAM YOU.
% SETS BASAL OSCILLATIONS, CAN OVERWRITE DOWNSTREAM
B(:,1) = BMaxODE * ones(N_r,1);
A(:,1) = AMaxODE * ones(N_r,1);
M(:,1) = zeros(N_r,1);

%% DEFINE ddr matrices; implicitly have Neumann programmed in

ddr = zeros(N_r,N_r);
for i=2:1:N_r-1
    ddr(i,i-1) = -1;
    ddr(i,i+1) = 1;
end
ddr = 1 / (2*dx) * ddr;

%% DEFINE d2dr2 matrix; implicitly have Neumann programmed in
d2dr2 = zeros(N_r,N_r);
for i=2:1:N_r-1
    d2dr2(i,i-1) = 1;
    d2dr2(i,i) = -2;
    d2dr2(i,i+1) = 1;
end
d2dr2 = 1 / (dx.^2) * d2dr2;



%% RUN WITH GRADIENT I.C.: DIFFUSION
%% SIMULATION RUN PARAMETERS
% step sizes
dr = 1;
% Declare an isotropic spacer, dx
dx = dr;
dt = 0.1;
% limits; small size for replicating figures
R_dim = 40; %microns
T_dim = 1000;
% N, +2 denotes phantom points to help with Neumann B.C.'s
N_r = round(R_dim / dr) + 2;
N_t = round(T_dim / dt);
% Dimensionalize time
tvec = 0:1:N_t-1;
tvec = 5 * dt * tvec;
%% CONSTANTS
epsilon = 0.1;
K = 1;
eta_B = 1;
% eta_M and eta_A not clearly labeled
eta_M = 1;
eta_A = 1;
E = 0.1;
B0_c = 4;
theta = 0.2;
noise = 0.01;
% Additional diffusion term for VASP; define based on epsilon2 =
% D_A/cell_l^2/gamma
% VASP estimated diffusivity
D_A = 0.2; %um/sec
cell_l = 40; %um
gamma = 1; % sec
% Non-dimensionalized lateral diffusivity
epsilon2 = D_A / cell_l^2 / gamma * 2;
%% INITIAL CONDITIONS; refer to icy.m for usage
BMaxODE = 8.80342361002194;
BMinODE = 2.46024797933793;
AMaxODE = 5.79255834170108;
AMinODE = 2.26789894792708;
MMaxODE = 2.69213920749395;
MMinODE = 0.677418361848779;
B = noise * rand(N_r,N_t);
A = zeros(N_r,N_t);
M = zeros(N_r,N_t);
V = zeros(N_r,N_t);

% INITIAL CONDITIONS SENT BY AUTHORS
% DON'T YOU DARE TOUCH THESE I.C.s FOR BAM OR I'LL BAM YOU.
% SETS BASAL OSCILLATIONS, CAN OVERWRITE DOWNSTREAM
B(:,1) = BMaxODE * ones(N_r,1);
A(:,1) = AMaxODE * ones(N_r,1);
M(:,1) = zeros(N_r,1);
% RUN WITH BOUNDARY CONDITIONS; targeted plots
delta = 0.8;
R = 0.12;
B = noise * rand(N_r,N_t);
A = zeros(N_r,N_t);
M = zeros(N_r,N_t);
V = zeros(N_r,N_t);
% INITIAL CONDITIONS FOR GRADIENT
B(:,1) = BMaxODE * ones(N_r,1);
A(:,1) = AMaxODE * ones(N_r,1);
M(:,1) = zeros(N_r,1);
for i = round(N_r * 0.25):1:round(N_r * 0.75)
   B(i,1) = BMinODE;
   A(i,1) = AMinODE; 
   M(i,1) = MMaxODE;
end
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
    noise_B = noise * 2 * (rand(N_r,1) - ones(N_r,1));
    B(:,t+1) = B(:,t) + dt/epsilon * (diffusion_B + accumulation_B + loss_B + noise_B);
    
    %% CALCULATE A
    accumulation_A = ones(N_r,1) * delta;
    % B.C. Neumann
    A(1,t) = A(3,t);
    A(N_r,t) = A(N_r - 2,t);
    diffusion_A = epsilon2 * d2dr2 * A(:,t);
    loss_A = -1 ./ (ones(N_r,1) + M(:,t) + K .* B(:,t)) .* (ones(N_r,1) + eta_A .* V(:,t) + eta_M .* M(:,t) .* V(:,t)) .* A(:,t);
    A(:,t+1) = A(:,t) + dt * (accumulation_A + loss_A + diffusion_A);
    %% CALCULATE M
    accumulation_M = R * B(:,t);
    loss_M = -1 * (ones(N_r,1) * theta + eta_M * V(:,t)) .* M(:,t);
    M(:,t+1) = M(:,t) + dt * (accumulation_M + loss_M);
end
V_avg = mean(mean(V));
V_avg_plot(R_ct,delta_ct) = V_avg;
V_avg_plot_d(R_ct,delta_ct) = delta;
V_avg_plot_R(R_ct,delta_ct) = R;
delta_ct = delta_ct + 1;
V_norm = V/max(max(V));
B_norm = B/max(max(B));
A_norm = A/max(max(A));
M_norm = M/max(max(M));
rad_coor = round(size(V_norm,1)/2);

figure
imagesc(V(:,1:N_t-1))
xlabel('Time (0.5 Seconds)');
ylabel('Radial Position (10^{th} of a Micron)');
title('Waving Behavior Velocity Kymograph /w VASP diffusion');
colormap('jet')
colorbar
figure
imagesc(A(:,1:N_t-1))
xlabel('Time (0.5 Seconds)');
ylabel('Radial Position (10^{th} of a Micron)');
title('Waving Behavior VASP Kymograph /w VASP diffusion');
colormap('jet')
colorbar
V_diffusion = V;

%% TOGGLE ON/OFF TO MAKE VIDEOS (TAKES AWHILE TO RUN)


%% TRACK DISTANCE AFTER 1st WAVE (ARTIFICIAL EFFECTS DUE TO I.C.'s)
dista = zeros(size(V,1),N_t);
cell(:,1)=zeros(size(dista,1),1);
cell(:,2)=[1:size(dista,1)]';
toffset = 300;
for ii=[toffset:1:N_t]
    if ii > 1
        dista(:,ii) = dista(:,ii-1) + V_norm(:,ii-1) * dt;
    end
end
% Scale to micron
dista = 0.2 * dista;
figure
video1 = VideoWriter('lamellopodia_with_diffusion','MPEG-4');
open(video1);
for frame=1:25:size(dista,2)
%     cell(:,1)=cell(:,1)-sum(dista(:,frame:frame+24),2)/(dt*25);
%     plot(cell(:,1),cell(:,2),'k-');
    plot(dista(:,frame),cell(:,2),'k-')
    title('Leading Edge of Cell Simulated Movement');
    ylabel('Distance along leading Edge (microns)')
    xlabel('Position (microns)')
    xlim([0 100]);
    ylim([0 N_r]);
    vidframe = getframe(gcf);
    writeVideo(video1,vidframe);
end
close(video1);
%% RUN WITH GRADIENT I.C.: NO DIFFUSION
% RUN WITH BOUNDARY CONDITIONS; targeted plots
delta = 0.8;
R = 0.12;
B = noise * rand(N_r,N_t);
A = zeros(N_r,N_t);
M = zeros(N_r,N_t);
V = zeros(N_r,N_t);
% INITIAL CONDITIONS FOR GRADIENT
B(:,1) = BMaxODE * ones(N_r,1);
A(:,1) = AMaxODE * ones(N_r,1);
M(:,1) = zeros(N_r,1);
for i = round(N_r * 0.25):1:round(N_r * 0.75)
   B(i,1) = BMinODE;
   A(i,1) = AMinODE; 
   M(i,1) = MMaxODE;
end
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
    noise_B = noise * 2 * (rand(N_r,1) - ones(N_r,1));
    B(:,t+1) = B(:,t) + dt/epsilon * (diffusion_B + accumulation_B + loss_B + noise_B);
    
    %% CALCULATE A
    accumulation_A = ones(N_r,1) * delta;
    % B.C. Neumann
    A(1,t) = A(3,t);
    A(N_r,t) = A(N_r - 2,t);
    diffusion_A = 0 * d2dr2 * A(:,t);
    loss_A = -1 ./ (ones(N_r,1) + M(:,t) + K .* B(:,t)) .* (ones(N_r,1) + eta_A .* V(:,t) + eta_M .* M(:,t) .* V(:,t)) .* A(:,t);
    A(:,t+1) = A(:,t) + dt * (accumulation_A + loss_A + diffusion_A);
    %% CALCULATE M
    accumulation_M = R * B(:,t);
    loss_M = -1 * (ones(N_r,1) * theta + eta_M * V(:,t)) .* M(:,t);
    M(:,t+1) = M(:,t) + dt * (accumulation_M + loss_M);
end
%% TOGGLE ON/OFF TO MAKE VIDEOS (TAKES AWHILE TO RUN)
%% TRACK DISTANCE AFTER 1st WAVE (ARTIFICIAL EFFECTS DUE TO I.C.'s)
dista = zeros(size(V,1),N_t);
toffset = 300;
for ii=[toffset:1:N_t]
    if ii > 1
        dista(:,ii) = dista(:,ii-1) + V_norm(:,ii-1) * dt;
    end
end
% Scale to micron
dista = 0.2 * dista;
figure
video1 = VideoWriter('lamellopodia_no_diffusion','MPEG-4');
open(video1);
for frame=1:25:size(dista,2)
%     cell(:,1)=cell(:,1)-sum(dista(:,frame:frame+24),2)/(dt*25);
%     plot(cell(:,1),cell(:,2),'k-');
    plot(dista(:,frame),cell(:,2),'k-')
    title('Leading Edge of Cell Simulated Movement');
    ylabel('Distance along leading Edge (microns)')
    xlabel('Position (microns)')
    xlim([0 100]);
    ylim([0 N_r]);
    vidframe = getframe(gcf);
    writeVideo(video1,vidframe);
end
close(video1);
V_avg = mean(mean(V));
V_avg_plot(R_ct,delta_ct) = V_avg;
V_avg_plot_d(R_ct,delta_ct) = delta;
V_avg_plot_R(R_ct,delta_ct) = R;
delta_ct = delta_ct + 1;
V_norm = V/max(max(V));
B_norm = B/max(max(B));
A_norm = A/max(max(A));
M_norm = M/max(max(M));

rad_coor = round(size(V_norm,1)/2);
figure
imagesc(V(:,1:N_t-1))
xlabel('Time (0.5 Seconds)');
ylabel('Radial Position (10^{th} of a Micron)');
title('Waving Behavior Velocity Kymograph without VASP diffusion');
colormap('jet')
colorbar

figure
imagesc(A(:,1:N_t-1))
xlabel('Time (0.5 Seconds)');
ylabel('Radial Position (10^{th} of a Micron)');
title('Waving Behavior VASP Kymograph without VASP diffusion');
colormap('jet')
colorbar
V_orig = V;

%% Changes over time
figure
V_difference = V_orig - V_diffusion;
imagesc(V_difference(:,1:N_t-1))
xlabel('Time (0.5 Seconds)');
ylabel('Radial Position (10^{th} of a Micron)');
title('Velocity_d_i_f_f_u_s_i_o_n - Velocity_o_r_i_g_i_n_a_l');
colormap('jet')
colorbar

%% PLOT THE RMS
for tstep=1:1:N_t
    divisor = 1 / (N_r - 2);
    RMS(tstep) = divisor * sqrt(sum(V_difference(2:N_r-1,tstep).^2));
end
% Since Velocity Max = 1, scale to %
RMS = RMS * 100;
figure
plot(tvec,RMS)
title('Velocity RMS error due to Diffusion')
ylabel('Velocity RMS error Percent of Max')
xlabel('Time (seconds)')
grid on
ylim([0 100])