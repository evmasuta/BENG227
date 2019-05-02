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
D_A = 0.2; %um/sec
cell_l = 40; %um
gamma = 1; % sec
% Non-dimensionalized lateral diffusivity
epsilon2 = D_A / cell_l^2 / gamma;
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
    d2dr2(i,i-1) = -1;
    d2dr2(i,i) = 2;
    d2dr2(i,i+1) = -1;
end
d2dr2 = 1 / (dx.^2) * d2dr2;

delta_ct = 1;
R_ct = 1;
R_iter = 0:0.02:0.24;
d_iter = 0.6:0.05:1;
%% RUN WITH BOUNDARY CONDITIONS; targeted plots
for R = R_iter
    delta_ct = 1;
    for delta = d_iter
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
        tvec = 0:1:N_t-2;
        rad_coor = round(size(V_norm,1)/2);
        % Dimensionalize time
        tvec = 5 * tvec * dt;
        if R == 0.24 && delta == 0.6
            figure
            plot(tvec,V_norm(rad_coor,1:N_t-1),tvec,B_norm(rad_coor,1:N_t-1),...
                tvec,A_norm(rad_coor,1:N_t-1),tvec,M_norm(rad_coor,1:N_t-1),'LineWidth',3)
            legend('V','B','A','M');
            xlabel('Time (Seconds)');
            ylabel('Normalized Quantity');
            ylim([0 1.1])
            grid on
            title(['Stalled Behavior at a Radial Coordinate of ',...
                num2str(rad_coor/dr),' Microns']);
            figure
            imagesc(V(:,1:N_t-1))
            xlabel('Time (Seconds)');
            ylabel('Radial Position (10^{th} of a Micron)');
            title('Stalled Behavior Velocity Kymograph');
            colormap('jet')
            colorbar
        elseif R == 0 && delta == 1
            figure
            plot(tvec,V_norm(rad_coor,1:N_t-1),tvec,B_norm(rad_coor,1:N_t-1),...
                tvec,A_norm(rad_coor,1:N_t-1),tvec,M_norm(rad_coor,1:N_t-1),'LineWidth',3)
            legend('V','B','A','M');
            xlabel('Time (10^{th} of a Second)');
            ylabel('Normalized Quantity');
            ylim([0 1.1])
            grid on
            title(['Smoothly Motile Behavior at a Radial Coordinate of ',...
                num2str(rad_coor/dr),' Microns']);
            figure
            imagesc(V(:,1:N_t-1))
            xlabel('Time (Seconds)');
            ylabel('Radial Position (10^{th} of a Micron)');
            title('Smoothly Motile Behavior Velocity Kymograph');
            colormap('jet')
            colorbar
        elseif R == 0.12 && delta == 0.8
            figure
            plot(tvec,V_norm(rad_coor,1:N_t-1),tvec,B_norm(rad_coor,1:N_t-1),...
                tvec,A_norm(rad_coor,1:N_t-1),tvec,M_norm(rad_coor,1:N_t-1),'LineWidth',3)
            legend('V','B','A','M');
            xlabel('Time (10^{th} of a Second)');
            ylabel('Normalized Quantity');
            ylim([0 1.1])
            grid on
            title(['Waving Behavior at a Radial Coordinate of ',...
                num2str(rad_coor/dr),' Microns']);
            figure
            imagesc(V(:,1:N_t-1))
            xlabel('Time (Seconds)');
            ylabel('Radial Position (10^{th} of a Micron)');
            title('Waving Behavior Velocity Kymograph');
            colormap('jet')
            colorbar
        end
    end
    R_ct = R_ct+1;
    
end

%% Plot heatmap resembling figure 6
% Flip plot up down to pictographically mirror the figure
V_avg_plot = flipud(V_avg_plot);
V_avg_plot_d = flipud(V_avg_plot_d);
V_avg_plot_R = flipud(V_avg_plot_R);
figure()
imagesc(V_avg_plot)
colormap('jet')
colorbar
ax = gca;
ax.YTick = 1:1:length(R_iter);
ax.YTickLabel = fliplr(R_iter);
ax.XTick = 1:1:length(d_iter);
ax.XTickLabel = d_iter;
xlabel('VASP delivery rate (scaled units)')
ylabel('Adh. maturation rate (scaled units)')
title('Average velocity (normalized)')








