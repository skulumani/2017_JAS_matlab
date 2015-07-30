function [L_points, gamma] = libration_points(mu)

% 16 september 2014 -  added output for gamma values
%                 L4
%
% -L3------M1--+-----L1--M2--L2-
%
%                 L5


mu1 = 1 - mu;

%        x^5    x^4       x^3      x^2   x    const
poly1 = [1   -1*(3-mu)  (3-2*mu)  -mu   2*mu  -mu ];
poly2 = [1      (3-mu)  (3-2*mu)  -mu  -2*mu  -mu ];
poly3 = [1      (2+mu)  (1+2*mu)  -mu1 -2*mu1 -mu1];

% solve for roots of quintic polynomial

rt1 = roots(poly1);
rt2 = roots(poly2);
rt3 = roots(poly3);



for k=1:5
        if isreal(rt1(k)) gam1=rt1(k); end
        if isreal(rt2(k)) gam2=rt2(k); end
        if isreal(rt3(k)) gam3=rt3(k); end
end

L1 = (1-mu) - gam1;
L2 = (1-mu) + gam2;
L3 = -mu - gam3;




L1 = [L1,0];
L2 = [L2,0];
L3 = [L3,0];
L4 = [1/2-mu,1/2*sqrt(3)];
L5 = [1/2-mu, -1/2*sqrt(3)];
L_points = [L1;L2;L3;L4;L5];
gamma = [gam1;gam2;gam3];