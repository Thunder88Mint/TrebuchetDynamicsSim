% This script is to be used to test your function to simulates the motion
% of a trebuchet. You are welcome to change model parameters, initial
% conditions, etc., but your function MUST work with an unmodified version
% of this script in order to be considered fully correct. Pay special
% attention to initial conditions and units.

clear;
clc;

% ---------------------------------------------------------------------
% Model Parameters
% ---------------------------------------------------------------------

% Parameters that are sent to your function
mp = 30;        % Projectile mass (kg)
mw = 3000;      % Counterweight mass (kg)
dw = 1.5;       % Distance to counterweight (m)

% Parameters used in determining initial conditions
dr = 8;         % Distance to tip of arm (m)
l = 6;          % Rope length (m)
h = 6;          % Height of pivot (m)

% ---------------------------------------------------------------------
% Initialization
% ---------------------------------------------------------------------

% Simulation time -- Your function must work with whatever time vector I
% send to it
% t = 0:0.0001:12;
t = 0:0.01:12;

% Initial conditions -- Note that specifying the initial angle of the arm
% uniquely defines the initial x
theta0 = -45;
phi0 = 0;
x0 = dr*cosd(theta0) - sqrt(l^2 - (dr*sind(theta0) + h)^2);
y0 = -h;

% ---------------------------------------------------------------------
% Simulation
% ---------------------------------------------------------------------

% Call your function -- Remember that your function and its file must be
% named with your last name, all lowercase. Everything you provide should
% be contained in a single file. There can be multiple functions within
% that file, but this driver function will only call the function called
% your_last_name().
[theta,phi,x,y,N,T] = greener(theta0,phi0,x0,y0,t,mp,mw,dw);

% ---------------------------------------------------------------------
% Plots
% ---------------------------------------------------------------------

% Motion vs. Time
figure(1);
subplot(4,1,1);
plot(t,x);
ylabel('x (m)');
subplot(4,1,2);
plot(t,y);
ylabel('y (m)');
subplot(4,1,3);
plot(t,theta);
ylabel('\theta (deg)');
subplot(4,1,4);
plot(t,phi);
ylabel('\phi (deg)');
xlabel('t (s)');

% Constraint Forces vs. Time
figure(2);
subplot(2,1,1);
plot(t,N);
ylabel('N (N)');
subplot(2,1,2);
plot(t,T);
ylabel('T (N)');
xlabel('t (s)');

% Projectile Motion (y vs. x)
figure(3);
plot(x,y);
ylabel('y (m)');
xlabel('x (m)');