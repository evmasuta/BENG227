%% BENG 227 MIDTERM PROJECT Evan Masutani
%% HOUSEKEEPING
clear all
close all
clc
%% INITIAL CONDITIONS
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3679584/
c_0 = 0.1*10^-6; %M
ce_0 = 50 * 10^-6; %M
p_0 = 0*10^-6; %M
h_0 = 0;
%% ODE45 SOLVER
tspan = [0 1000];
[t,y] = ode45(@fluxes,tspan,[c_0;ce_0;p_0;h_0]);
figure
subplot(2,2,1)
plot(t,y(:,1))
title('Cytosol Ca2+')
subplot(2,2,2)
plot(t,y(:,2))
title('ER Ca2+')
subplot(2,2,3)
plot(t,y(:,3))
title('IP3')
subplot(2,2,4)
plot(t,y(:,4))
title('h')