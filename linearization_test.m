% 3 April 2015 - compare linearized system to nonlinear equations of motion
clear all
clc
close all

constants = crtbp_constants;
constants.um = 0.005;

step_exp = 4;
num_steps = round(1*10^step_exp)+1;
% num_steps = 300+1;

[L_points, ~] = libration_points(constants.mu);
constants.center_vec = L_points(1,:);

constants.alpha_d = 0*pi/180 ; % angle that defines the poincare section
% constants.traj_fig = figure(1);
% generate an initial periodic orbit
% initial_energy = constants.e_desired;
% [x0_i, T_i, E_i, phi_cross_i] = periodic_orbit_pcrtbp(constants.l_point, constants.e_desired, constants);
x0_i = [0.815614054266804 0 0 0.192227407664904]';
T_i = 1.407478324303006;

% generate a target periodic orbit
% target_energy = initial_energy + 0.005;
% [x0_t, T_t, E_t, phi_cross_t] = periodic_orbit_pcrtbp(constants.l_point, target_energy, constants);
x0_t = [0.813217758536162 0 0 0.219019807293238];
T_t = 1.429672365611922;
% propogate the intial and target orbits to the poincare section using
% variational integrator

% generate initial guess of state and costate histories
h0 = [0.003042140861866;... % hr_x
    0.005574701108068;... % hr_y
    0.002230056745460;... % hv_x
    -0.001157375275785]; % hv_y

% x0 = rand(4,1);
x0 = x0_i;
% h0 = rand(4,1);

deltax = 1e-4*rand(size(x0));
deltah = 1e-4*rand(size(h0));

tf = 2;
t = linspace(0,tf,num_steps);
h = t(2)-t(1);
num_states = size(x0,1);

constants.control_switch = 'on';
[~, state_guess,~] = pcrtbp_optimal_var(x0,h0,0,tf, 'trap', num_steps, constants);
[~, state_delta,~] = pcrtbp_optimal_var(x0+deltax,h0+deltah,0,tf, 'trap', num_steps, constants);
% initial guess of state and costate histories
x_i = state_guess(:,1:4);
h_i = state_guess(:,5:8);

% calculate state history loop using state transition matrix
A = zeros(2*num_states,2*num_states,num_steps);
e = zeros(2*num_states,num_steps);
z = zeros(num_steps,2*num_states);

approx = zeros(2*num_states,num_steps);
% define the initial states
% [dfdx, dfdh, dgdx, dgdh, f, g] = pcrtbp_trap_lin(x_i(1,:),h_i(1,:),constants.mu,h); % linearize system about current trajectory
z(1,:) = [x0;h0]';
% A(:,:,1) = [dfdx dfdh;dgdx dgdh];
% e(:,1) = - A(:,:,1)*[x_i(1,:)';h_i(1,:)'] + [f;g];
for k=1:num_steps-1
    [dfdx, dfdh, dgdx, dgdh, f, g] = pcrtbp_trap_lin(x_i(k,:),h_i(k,:),h,constants); % linearize system about current trajectory
    
    A(:,:,k) = [dfdx dfdh;dgdx dgdh];

    delta_state = abs(state_delta(k,:)-state_guess(k,:));
    
    approx(:,k) = A(:,:,k)*delta_state';
    
    e(:,k) = - A(:,:,k)*[x_i(k,:)';h_i(k,:)'] + [f;g];

    z(k+1,:)=(A(:,:,k)*z(k,:)'+e(:,k))';

end

x = z(:,1:num_states);
lambda = z(:,num_states+1:end);

% solve for N
% calculate Psi
Psi = eye(size(A(:,:,1)));

for ii = 1:num_steps-1
    Psi = A(:,:,ii)*Psi;
end

% calculate particular solution
% get the product of the A matrices and store in another array
M = zeros(size(A));
M(:,:,1) = eye(size(A(:,:,1)));

for ii = 1:num_steps-1
    M(:,:,ii+1) = M(:,:,ii)*A(:,:,num_steps-ii);
end

% multiply M x e and add them up
P = zeros(size(e(:,1)));
for ii = 1:num_steps-1
    P = P+M(:,:,num_steps-ii)*e(:,ii);
end

zn = Psi*z(1,:)' + P;

zn-z(end,:)'
% compare perturbations
figure

subplot(4,2,1)
plot(t,abs(state_delta(:,1)-state_guess(:,1)),t,approx(1,:))
title('State Perturbation Comparison')
ylabel('x')
legend('f(x^*+\delta x) - f(x^*)','DF \delta')
grid on
subplot(4,2,3)
plot(t,abs(state_delta(:,2)-state_guess(:,2)),t,approx(2,:))
ylabel('y')
grid on
subplot(4,2,5)
plot(t,abs(state_delta(:,3)-state_guess(:,3)),t,approx(3,:))
ylabel('xd')
grid on
subplot(4,2,7)
plot(t,abs(state_delta(:,4)-state_guess(:,4)),t,approx(4,:))
ylabel('yd')
xlabel('t')
grid on

subplot(4,2,2)
plot(t,abs(state_delta(:,5)-state_guess(:,5)),t,approx(5,:))
title('Costate Pertubtion Difference')
ylabel('hx')
grid on
subplot(4,2,4)
plot(t,abs(state_delta(:,6)-state_guess(:,6)),t,approx(6,:))
ylabel('hy')
grid on
subplot(4,2,6)
plot(t,abs(state_delta(:,7)-state_guess(:,7)),t,approx(7,:))
ylabel('hxd')
grid on
subplot(4,2,8)
plot(t,abs(state_delta(:,8)-state_guess(:,8)),t,approx(8,:))
ylabel('hyd')
xlabel('t')
grid on


figure
grid on
hold all;
plot(state_guess(:,1),state_guess(:,2))
plot(x(:,1),x(:,2))
xlabel('x')
ylabel('y')
legend('NL EOM (1)','Lin EOM (4)')

% compare state trajectories
figure
grid on
hold on
subplot(4,2,1)
plot(t,state_guess(:,1),t,x(:,1))
title('State Comparison')
ylabel('x')
legend('NL EOM (1)','Lin EOM (4)')
grid on
subplot(4,2,3)
plot(t,state_guess(:,2),t,x(:,2))
ylabel('y')
grid on
subplot(4,2,5)
plot(t,state_guess(:,3),t,x(:,3))
ylabel('xd')
grid on
subplot(4,2,7)
plot(t,state_guess(:,4),t,x(:,4))
ylabel('yd')
xlabel('t')
grid on

subplot(4,2,2)
plot(t,abs(state_guess(:,1)-x(:,1)))
title('Difference')
grid on
subplot(4,2,4)
plot(t,abs(state_guess(:,2)-x(:,2)))
grid on
subplot(4,2,6)
plot(t,abs(state_guess(:,3)-x(:,3)))
grid on
subplot(4,2,8)
plot(t,abs(state_guess(:,4)-x(:,4)))
xlabel('t')
grid on

% compare costate trajecoties
figure
grid on
hold on
subplot(4,2,1)
plot(t,state_guess(:,5),t,lambda(:,1))
title('Costate comparison')
ylabel('x')
legend('NL EOM (1)','Lin EOM (4)')
grid on
subplot(4,2,3)
plot(t,state_guess(:,6),t,lambda(:,2))
ylabel('y')
grid on
subplot(4,2,5)
plot(t,state_guess(:,7),t,lambda(:,3))
ylabel('xd')
grid on
subplot(4,2,7)
plot(t,state_guess(:,8),t,lambda(:,4))
ylabel('yd')
xlabel('t')
grid on

subplot(4,2,2)
plot(t,abs(state_guess(:,5)-lambda(:,1)))
title('Difference')
grid on
subplot(4,2,4)
plot(t,abs(state_guess(:,6)-lambda(:,2)))
grid on
subplot(4,2,6)
plot(t,abs(state_guess(:,7)-lambda(:,3)))
grid on
subplot(4,2,8)
plot(t,abs(state_guess(:,8)-lambda(:,4)))
xlabel('t')
grid on
