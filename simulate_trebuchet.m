function [theta,phi,x,y,N,T] = simulate_trebuchet(theta0,phi0,x0,y0,t,mp,mw,dw)

% Unit Convert Inputs
theta0 = theta0 *pi/180;    % rad
phi0 = phi0 *pi/180;        % rad


% Model parameters & Constants
ma = 150;       % Arm Mass (kg)
I_AO = 1466;    % Moment of Inertia of arm about O (kg m2)
I_WG = 6750;    % Moment of Inertia of counterweight about center of mass (kg m2)
dr = 8;         % (m)
da = 1.5;       % (m)
dg = 1.5;       % (m)
mu = 2;         % Friction Coeff - Ground and projectile
g = 9.81;       % gravity (m/s2)


% Initial Conditions
q0 = [x0; y0; theta0; phi0; 0; 0; 0; 0];  


% ODE Solver
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15);
[t,z] = ode45(@(t,x) trebuchet(t,x,mp,ma,mw,I_AO,I_WG,dr,da,dw,dg,mu,g), t, q0,options);


% Call function again to extract lagrange multipliers
lambda1 = zeros(length(t),1);
lambda2 = zeros(length(t),1);
for i = 1:length(t)
    [~, lambda] = trebuchet(t(i),z(i,:),mp,ma,mw,I_AO,I_WG,dr,da,dw,dg,mu,g);
    lambda1(i) = lambda(1);
    lambda2(i) = lambda(2);
end


% Extract variables of interest
x = z(:,1);
y = z(:,2);
theta = z(:,3);
phi = z(:,4);
xdot = z(:,5);
ydot = z(:,6);
thdot = z(:,7);
phidot = z(:,8);


% Calculate generalized constraint forces
% Calculate Tension Force
T = -lambda1.*sqrt((xdot+dr*thdot.*sind(theta)).^2 + (ydot-dr*thdot.*cosd(theta)).^2);
% Calculate Normal Force
N = lambda2;


% Unit Convertions before output
theta = theta * 180/pi;
phi = phi * 180/pi;



%% My Plots - Troubleshooting

% figure(1)
% plot(x,y,'c',0,0,'o')
% hold on
% % Arm
% plot([0 dr*cosd(theta(end))], [0 dr*sind(theta(end))], 'm')
% % Rope
% % plot([x(end), dr*cosd(theta(end))], [y(end) dr*sind(theta(end))], 'k')
% hold off
% xlabel('x distance (m)')
% ylabel('y distance (m)')
% 
% % plot(t,x)
% 
% figure(2)
% subplot(4,1,1)
% plot(t,xdot)
% xlabel('Time (s)')
% ylabel('xdot (m/s)')
% subplot(4,1,2)
% plot(t,ydot)
% xlabel('Time (s)')
% ylabel('ydot (m/s)')
% subplot(4,1,3)
% plot(t,thdot)
% xlabel('Time (s)')
% ylabel('thetadot (rad/s)')
% subplot(4,1,4)
% plot(t,phidot)
% xlabel('Time (s)')
% ylabel('phidot (rad/s)')
% 
% 
% figure(3)
% plot(t(1:length(N)),N)
% ylabel('Normal Force (N)')
% 
% figure(4)
% plot(t(1:length(T)),T)
% ylabel('Tension Force (N)')
% 
% 
% % Plot Constraints
% trans1_ind = find(N <= 0, 1, 'first');
% figure(5)
% % trans1_ind = find(N == 0, 1);
% plot(t(1:trans1_ind), y(1:trans1_ind))
% ylabel('y-coordinate')
% title('Ground Constraint Check')
% 
% figure(6)
% trans2_ind = find(T(2:end) <= 0, 1, 'first');
% l = sqrt((x(1:trans2_ind)-dr*cosd(theta(1:trans2_ind))).^2 + (y(1:trans2_ind)-dr*sind(theta(1:trans2_ind))).^2);
% plot(t(1:trans2_ind), l)
% ylabel('Rope Length')
% title('Rope Constraint Check')

end


function [xdot,lambda] = trebuchet(t,x,mp,ma,mw,I_AO,I_WG,dr,da,dw,dg,mu,g)
    persistent N phase

    % Reset global variables
    if t == 0
        N = mp*g;
        phase = 1;
    end

    % Extract generalized coordinates
    q1 = x(1); % x
    q2 = x(2); % y
    q3 = x(3); % theta
    q4 = x(4); % phi
    q1dot = x(5); % xdot
    q2dot = x(6); % ydot
    q3dot = x(7); % theta_dot
    q4dot = x(8); % phi_dot

    % Mass Matrix [M]
    M = [mp 0 0 0;...
        0 mp 0 0;...
        0 0 I_AO+mw*dw^2 mw*dw*dg*sin(q3-q4);
        0 0 mw*dw*dg*sin(q3-q4) I_WG+mw*dg^2];

    % Force Matrix [F]
    F = [-mu*N;
        -mp*g;
        mw*dw*dg*q4dot^2*cos(q3-q4) + g*cos(q3)*(mw*dw - ma*da);
        -mw*dw*dg*q3dot^2*cos(q3-q4) - mw*g*dg*sin(q4)];

    % Vector of generalized velocities
    qdot = [q1dot; q2dot; q3dot; q4dot];


    % Constraint Matrix [a] & da/dt
    % Solvers
    switch phase
        case 1 % Sliding
            a = [q1-dr*cos(q3) q2-dr*sin(q3) dr*(q1*sin(q3)-q2*cos(q3)) 0;
                0 1 0 0];
            dadt = [q1dot+dr*q3dot*sin(q3) q2dot-dr*q3dot*cos(q3) dr*(q1dot*sin(q3)+q1*q3dot*cos(q3)-q2dot*cos(q3)+q2*q3dot*sin(q3)) 0;
                0 0 0 0];

            % Solver
            qddot_lambda = [M -a'; -a zeros(2)] \ [F; dadt*qdot];

            % Extract generalized accelerations and lagrange multipliers
            qddot = qddot_lambda(1:4);
            lambda = qddot_lambda(5:end);
        case 2 % Swinging
            a = [q1-dr*cos(q3) q2-dr*sin(q3) dr*(q1*sin(q3)-q2*cos(q3)) 0];
            dadt = [q1dot+dr*q3dot*sin(q3) q2dot-dr*q3dot*cos(q3) dr*(q1dot*sin(q3)+q1*q3dot*cos(q3)-q2dot*cos(q3)+q2*q3dot*sin(q3)) 0];

            % Solver
            qddot_lambda = [M -a'; -a zeros(1)] \ [F; dadt*qdot];

            % Extract generalized accelerations and lagrange multipliers
            qddot = qddot_lambda(1:4);
            lambda = [qddot_lambda(5:end); 0];
        case 3 % Flying
            qddot = M\F;
            lambda = [0; 0];
    end


    % Return state derivatives (Output)
    xdot = [qdot; qddot];


    % Evaluate Phase Transitions
    switch phase
        case 1 % Sliding to swinging
            N = lambda(2);
            if N <= 0
                phase = 2;
                N = 0;
            end

        case 2 % Swinging to Flying
            % calculate slopes
            k_arm = tan(q3);
            k_rope = (q2-dr*sin(q3)) / (q1-dr*cos(q3));
            if abs(k_arm - k_rope)/abs(k_arm) < 0.7
                phase = 3;
            end
    end
    

end





