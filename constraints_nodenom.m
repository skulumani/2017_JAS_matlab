% 19 June 2015
% modified constraint function that avoids the inverse tangent function and
% small values in the denominator
function [m, dmdx,dm1xdx,dm2xdx] = constraints_nodenom(statef,costatef,constants)
num_con = 2;
xf = statef(1);
yf = statef(2);
xdf = statef(3);
ydf = statef(4);

xcf = constants.xcf;
xcen = constants.center_vec;

alpha = constants.alpha_d;
theta = constants.theta_d;

m = zeros(num_con,1);

m(1) = sin(alpha)^2*((xcen(1) - xf)^2+(xcen(2) - yf)^2) - (xcen(2) - yf)^2;
m(2) = sin(theta)^2*((xcf(3) - xdf)^2+(xcf(1)-xf)^2) - (xcf(3)-xdf)^2;
    

dmdx = [-sin(alpha)^2*(2*xcen(1) - 2*xf), -2*(sin(alpha)^2 - 1)*(xcen(2) - yf),0, 0;...
        -sin(theta)^2*(2*xcf(1) - 2*xf),0, -2*(sin(theta)^2 - 1)*(xcf(3) - xdf), 0];
    
dm1xdx = ...
[ 2*sin(alpha)^2,                  0, 0, 0;...
              0, 2*sin(alpha)^2 - 2, 0, 0;...
              0,                  0, 0, 0;...
              0,                  0, 0, 0];


dm2xdx = ...
    [ 2*sin(theta)^2, 0,                  0, 0;...
        0, 0,                  0, 0;...
        0, 0, 2*sin(theta)^2 - 2, 0;...
         0, 0,                  0, 0];

end