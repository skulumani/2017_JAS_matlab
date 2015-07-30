function state_dot = pcrtbp_ode_update(state,mu,u)


% define some parameters

x = state(1);
y = state(2);

xd = state(3);
yd = state(4);

A = 2 * [0 1;...
    -1 0 ];

% constants/parameters

r1 = sqrt((x+mu)^2 + y^2 );
r2 = sqrt( (x+mu-1)^2+ y^2 );

ux = x - ((1-mu)*(x+mu))/(r1^3) - mu*(x+mu-1)/(r2^3);
uy = y - (1-mu)*y/r1^3  - mu *y/r2^3;

r_dot = [xd;yd];
v_dot = A*[xd;yd] + [ux;uy] + u;

state_dot = [r_dot;v_dot];