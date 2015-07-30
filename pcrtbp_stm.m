function state_dot = pcrtbp_stm(t,state,mu)

%propogate the EOMs and STM together for Planar CRTBP

x = state(1);
y = state(2);

xd = state(3);
yd = state(4);

phi = state(5:20);
phi = reshape(phi,4,4);

[Df_mat] = linearized_eom_mat([x;y;xd;yd], mu); % only need the 4x4 matrix

% Df_mat(6,:) = []; Df_mat(:,6) = [];
% Df_mat(3,:) =[]; Df_mat(:,3) = [];
 % upper 4x4 for planar case
 
 
 
% constants/parameters

d1 = sqrt((x+mu)^2 + y^2 );
d2 = sqrt( (x+mu-1)^2+ y^2 );

ux = x - ((1-mu)*(x+mu))/(d1^3) - mu*(x+mu-1)/(d2^3);
uy = y - (1-mu)/d1^3 * y - mu *y/d2^3;
% uz = -(1-mu)/d1^3 *z - mu*z/d2^3;

x_dot = [xd;yd;...
    2*yd + ux;...
    -2*xd + uy];

phidot = Df_mat * phi; % variational equation
phidot = reshape(phidot,16,1);

state_dot = [x_dot;phidot];