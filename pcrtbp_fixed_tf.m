% 16 dec 14 - modfying for use as optimal control code to a fixed final
% state

function [] = pcrtbp_fixed_tf(xc0, xcf)

constants.mu = 0.0125;
um = 0.75;
constants.plot_scale = 'nondim';
constants.l_scale = 384400; % distance of M2 from M1 (SMA) km

ode_options = odeset('RelTol',1e-9,'AbsTol',1e-9);
constants.ode_options = ode_options;

optfsolve = optimoptions(@fsolve,'Display','iter','TolFun',1e-9,'TolX',1e-9);
% tspan = [0 0.8];
% tspan = [0 1.444];
tspan = [0 1.3075]; % L1 transfer final
% tspan = [0 0.355837553671024]; % geo transfer final

% xc0 = [0.8352;0;0;0];
% xc0 = [-0.4127;0.1691;0.3565;-1.2079];
% xc0 = [0.8119;0;0;0.2342];
xc0 = [0.815614054266804 0 0 0.192227407664904]'; % L1 reach transfer
% xc0 = [0.17541888429434552000; 0.00000000000000000003; -0.40753131851731217000; 1.71883588181152410000 ]; % geo final transfer

constants.xc0 = xc0;

constants.um = um;
constants.tspan = tspan;
constants.ode_options = ode_options;
rand_theta = 0 + (2*pi-0).*rand;
constants.sub_opt = [cos(rand_theta); sin(rand_theta)];
% propogate xc0 to find xcf for optimization cost
% constants.xcf = [0.8352;0;0;0]';
% constants.xcf = [0.2;-0.3;0;0]';
constants.xcf = [0.927258145099694 0  -0.036751604638168 -0.369096025781395]; % L1 reach transfer final
% constants.xcf = [ 0.175319307882103 -0.000000000000020 -0.282163264918425 2.717676740320596]; % geo transfer final

hguess = [-0.0049;-0.0077;0.0007;-0.0060];
[h0,fval,exitflag,output] = fsolve(@(h0)obj(h0,tspan,constants),hguess,optfsolve);

constants.control_switch = 'on';
[t,state] = ode45(@(t,state)pcrtbp_ode(t,state,constants),constants.tspan,[xc0;h0],ode_options);

pos = state(:,1:2);
vel = state(:,3:4);
hr = state(:,5:6);
hv = state(:,7:8);
norm_hv = sqrt(sum(hv.^2,2));
u = -constants.um*(hv./repmat(norm_hv,1,2))';
u = u';

figure(1)
hold all
plot(pos(:,1),pos(:,2),'g','linewidth',4)
% plot_trajectories(t, [pos vel], energyconst(constants.xc0',constants.mu), figure(1), constants)
plot(constants.xc0(1),constants.xc0(2),'ko')
plot(constants.xcf(1),constants.xcf(2),'ko')

% quiver(pos(1:5:end,1),pos(1:5:end,2),u(1:5:end,1),u(1:5:end,2));
J=1/2*(state(end,1:4)' - constants.xcf(1:4)')'*[eye(2,2) zeros(2,2);zeros(2,2) zeros(2,2)]*(state(end,1:4)' - constants.xcf(1:4)')

figure 
grid on
hold on
plot(t, u(:,1), t, u(:,2))
% plot(t,hv(:,1),t, hv(:,2))
end

function state_dot=pcrtbp_ode(t,state,constants)
% 6 november 2014 - added control switch for reachability purposes

% seperate out the states

mu = constants.mu;
% define some parameters
mu1 = 1-mu; % larger primary
mu2 = mu; % smaller primary

x = state(1);
y = state(2);
% z = state(3);

xd = state(3);
yd = state(4);
% zd = state(6);

hx = state(5);
hy = state(6);
% hz = state(9);
hxd = state(7);
hyd = state(8);
% hzd = state(12);

hr = [hx;hy];
hv = [hxd;hyd];

% constants/parameters

r1 = sqrt((x+mu)^2 + y^2 );
r2 = sqrt( (x+mu-1)^2+ y^2 );

ux = x - ((1-mu)*(x+mu))/(r1^3) - mu*(x+mu-1)/(r2^3);
uy = y - (1-mu)*y/r1^3  - mu *y/r2^3;
% uz = -(1-mu)*z/d1^3  - mu*z/d2^3;

% compute the second partial derivatives of the gravitational potential
Uxx = 1 - mu1 * ( 1/r1^3 - ( 3*((x + mu2)^2)/r1^5 ) ) - mu2 * (1/r2^3 - 3*(x - mu1)^2 /r2^5);
Uyy = 1 - mu1 * ( 1/r1^3 -  3*y^2/r1^5 ) - mu2 * ( 1/r2^3 -3*y^2/r2^5);

Uxy = 3*mu1*y*(x+mu2)/r1^5 + 3*mu2*y*(x-mu1)/r2^5; Uyx = Uxy;


G = [Uxx Uxy ;...
    Uyx Uyy ];
A = 2 * [0 1;...
    -1 0 ];
B = [zeros(2,2); eye(2,2)];

switch constants.control_switch
    case 'off'
        u = zeros(2,1);
    case 'on'
        % express control vector in terms of costates
        u = -constants.um *hv/1;
%         norm_hv = norm(hv);
%         if norm_hv < 1e-9
%             u = zeros(2,1);
%         else
%              cs = sign(B'*[hr;hv]);
%             u = [0;0];
%             switch cs(1)
%                 case 1
%                     u(1) = -constants.um;
%                 case -1
%                     u(1) = constants.um;
%             end
%             
%             switch cs(2)
%                 case 1
%                     u(1) = -constants.um;
%                 case -1
%                     u(1) = constants.um;
%             end
%             u = -constants.um *hv/norm(hv);
% %                 u = -constants.um *sign(hv)/norm(hv);
%         end
    case 'sub'
        rand_vec = constants.sub_opt;
        u = constants.um*rand_vec/norm(rand_vec);
end


% u = -constants.acc_max *hv;
% derivatives
r_dot = [xd;yd];
v_dot = A*[xd;yd] + [ux;uy] + u;

hr_dot = - G'*hv;
hv_dot = -hr - A'*hv;

state_dot = [r_dot;v_dot;hr_dot;hv_dot];

end

function F = obj(h0,tspan,constants)

initial_state = [constants.xc0;h0];
% integrate the equations of motion for the given h_0 
constants.control_switch = 'on';
[~,state] = ode45(@(t,state)pcrtbp_ode(t,state,constants),tspan,initial_state,constants.ode_options);
% subtract the final state from the desired final state 
F = ( (state(end,1:4)'- constants.xcf') );


end
