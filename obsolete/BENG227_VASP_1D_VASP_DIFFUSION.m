%% BENG 227 MIDTERM PROJECT Evan Masutani
%% HOUSEKEEPING
clear all;
close all;
clc;
%% SIMULATION RUN PARAMETERS
% step sizes
dr = 0.1;
% Declare an isotropic spacer, dx
dx = dr;
dt = 0.01;
% limits
R_dim = 60; %microns
T_dim = 100;
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

D_A = 0.01;
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


% B(1:N_IC,1) = BMinODE * ones(N_IC,1);
% A(1:N_IC,1) = AMinODE * ones(N_IC,1);
% M(1:N_IC,1) = MMaxODE * ones(N_IC,1);


% % Gradient I.C.


% B(:,1) = BMinODE * ones(N_r,1);
% A(:,1) = AMinODE * ones(N_r,1);
% M(:,1) = MMaxODE * ones(N_r,1);
for i = 200:1:400
   B(i,1) = BMinODE;
   A(i,1) = AMinODE; 
   M(i,1) = MMaxODE;

end




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
    noise_B = noise * 2 * (rand(N_r,1) - ones(N_r,1));
    B(:,t+1) = B(:,t) + dt/epsilon * (diffusion_B + accumulation_B + loss_B + noise_B);

    %% CALCULATE A
    accumulation_A = ones(N_r,1) * delta;
    % B.C. Neumann
    A(1,t) = A(3,t);
    A(N_r,t) = A(N_r - 2,t);
    diffusion_A = D_A * d2dr2 * A(:,t);
    loss_A = -1 ./ (ones(N_r,1) + M(:,t) + K .* B(:,t)) .* (ones(N_r,1) + eta_A .* V(:,t) + eta_M .* M(:,t) .* V(:,t)) .* A(:,t);
    A(:,t+1) = A(:,t) + dt * (accumulation_A + loss_A + diffusion_A);
    %% CALCULATE M
    accumulation_M = R * B(:,t);
    loss_M = -1 * (ones(N_r,1) * theta + eta_M * V(:,t)) .* M(:,t);
    M(:,t+1) = M(:,t) + dt * (accumulation_M + loss_M);
end

%% Normalize Results
V_norm = V/max(max(V));
B_norm = B/max(max(B));
A_norm = A/max(max(A));
M_norm = M/max(max(M));
rad_coor = round(size(V_norm,1)/2);

%% Plot Results
% figure
% plot(tvec,M(10,:))
% figure
% plot(tvec,A(10,:))
% figure
% plot(tvec,B(10,:))


%% TRACK DISTANCE AFTER 1st WAVE (ARTIFICIAL EFFECTS DUE TO I.C.'s)
dista = zeros(size(V,1),N_t);
toffset = 2000;
for ii=[toffset:1:N_t]
    if ii > 1
        dista(:,ii) = dista(:,ii-1) + V_norm(:,ii-1) * dt;
    end
end
cell(:,1)=zeros(size(dista,1),1);
cell(:,2)=[1:size(dista,1)]';

%% TOGGLE ON/OFF TO MAKE VIDEOS (TAKES AWHILE TO RUN)
figure
video1 = VideoWriter('lamellopodia','MPEG-4');
open(video1);
for frame=1:25:size(dista,2)
%     cell(:,1)=cell(:,1)-sum(dista(:,frame:frame+24),2)/(dt*25);
%     plot(cell(:,1),cell(:,2),'k-');
    plot(dista(:,frame),cell(:,2),'k-')
    title('Leading Edge of Cell Simulated Movement');
    xlim([0 100]);
    ylim([0 N_r]);
    vidframe = getframe(gcf);
    writeVideo(video1,vidframe);
end
close(video1);
figure
imagesc(V(:,1:N_t-1))
xlabel('Time (10^{th} of a Second)');
ylabel('Radial Position (10^{th} of a Micron)');
title('Non-Dimensionalized Normalized Velocity Kymograph');
colorbar
figure
imagesc(B(:,1:N_t-1))
xlabel('Time (10^{th} of a Second)');
ylabel('Radial Position (10^{th} of a Micron)');
title('Non-Dimensionalized Barbed End Density Kymograph');
colorbar
% figure
% plot(tvec,V(10,:));
figure
imagesc(A(:,1:N_t-1))
xlabel('Time (10^{th} of a Second)');
ylabel('Radial Position (10^{th} of a Micron)');
title('Non-Dimensionalized VASP Concentration Kymograph');
colorbar
figure
imagesc(M(:,1:N_t-1))
xlabel('Time (10^{th} of a Second)');
ylabel('Radial Position (10^{th} of a Micron)');
title('Non-Dimensionalized Mature Adhesions Kymograph');
colorbar
figure
plot(tvec,V_norm(rad_coor,:),tvec,B_norm(rad_coor,:),...
    tvec,A_norm(rad_coor,:),tvec,M_norm(rad_coor,:),'LineWidth',3)
legend('V','B','A','M');
xlabel('Time (10^{th} of a Second)');
ylabel('Normalized Quantity');
title(['Non-Dimensionalized Quantities at a Radial Coordinate of ',...
    num2str(rad_coor/10),' Microns']);
% figure
% V_scale = V ./ max(max(V));
% imshow(V)