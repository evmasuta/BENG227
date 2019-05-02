%% BENG 227 MIDTERM PROJECT DIFFERENTIAL EQUATIONS
function dydt = fluxes(t,y)
%% BIOLOGICAL CONSTANTS
d = 1.5;
Kt = 0.1*10^-6; %M
Kc = 0.2*10^-6; %M
kf = 10; %Hz
Vserca = 0.9*10^-6; %M/s
g = 5.5;
Vpm = 0.11*10^-6; %M/s
a0 = 0.0027*10^-6; %M/s
Kp = 0.2*10^-6; %M
tmax = 1000; %Hz
tp = 0.027; %Hz
Kh = 0.08*10^-6; %M
Kbar = 1.9*10^-5;
Kserca = 0.2*10^-6; %M
Kce = 8*10^-6; %M
Kpm = 0.3*10^-6; %M
a1 = 0.07*10^-6; %M/s

%% PAPER-DERIVED VALUES
p_s = 0.12*10^-6; %M
kb = 0.4;

%% DIFFERENTIAL EQUATIONS
% Order is y(1) = c, y(2) = c_e, y(3) = p, y(4) = h

% FLUX TERMS
Jserca = Vserca*(y(1)^2 - Kbar*y(2)^2) / (y(1)^2 + Kserca^2);
Jin = a0 + a1*(Kce^4 / (Kce^4 + y(2)^4));
Jpm = Vpm*y(1)^2 / (Kpm^2 + y(1)^2);

% Jipr terms
A = 1 - y(3)^2 / (Kp^2 + y(3)^2);
B = 1 - A;
mbara = y(1)^4 / (Kc^4 + y(1)^4);
mbarb = mbara;
hbara = Kh^4 / (Kh^4 + y(1)^4);
hinf = hbara;
alpha = A*(1 - mbara*hbara);
beta = B*mbarb*y(4);
P0 = beta / (beta + kb*(beta + alpha));
Jipr = kf*P0*(y(2) - y(1));

% dcdt terms
dcdt = Jipr - Jserca + d*(Jin - Jpm);

% dc_edt terms
dc_edt = g*(Jserca - Jipr);

% dpdt terms
dpdt = tp*(p_s - y(3));

% dhdt terms
th = tmax * Kt^4 / (Kt^4 + y(1)^4);
dhdt = (hinf - y(4)) / th;

%% MODEL IP3 PULSE
if t > 356 && t < 357
   dpdt = dpdt + 2*p_s;
end
%% OUTPUTS

dydt = [dcdt; dc_edt; dpdt; dhdt];

end