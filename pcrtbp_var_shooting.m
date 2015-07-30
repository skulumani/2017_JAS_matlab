% 8 July 15 
% function to propogate the state, costate, and sensitivities at the same
% time


function [t, state,sens,H] = pcrtbp_var_shooting(x0,h0,Pxh0, Phh0, Pxx0, Phx0, t0,tf, mode, num_steps, constants)
t = linspace(t0,tf,num_steps)';

h = t(2)-t(1); % time step
mu = constants.mu;

% intialize arrays for the propogated states
state = zeros(num_steps,4);
costate = zeros(num_steps,4);

H = zeros(num_steps,1);

state(1,:) = x0';
costate(1,:) = h0';

switch mode
    case 'rect'
        
        % copy over the subfunction into here to avoid having multiple
        % subfunctions
        x = zeros(num_steps,1);
        y = zeros(num_steps,1);
        xd = zeros(num_steps,1);
        yd = zeros(num_steps,1);
        
        x(1) = x0(1);
        y(1) = x0(2);
        xd(1) = x0(3);
        yd(1) = x0(4);
        
        
        for k = 1:num_steps-1
            switch constants.control_switch
                case 'off'
                    u = zeros(2,1);
                case 'on'
                    % express control vector in terms of costates
                    hv = costate(k,3:4);
                    u = -constants.um *hv/norm(hv);
                case 'sub'
                    rand_vec = constants.sub_opt;
                    u = constants.um*rand_vec/norm(rand_vec);
            end
            
            h2 = h*h;
            h3 = h*h*h;
            
            r1 = sqrt((x(k)+mu)^2 + y(k)^2 );
            r2 = sqrt( (x(k)+mu-1)^2+ y(k)^2 );
            
            r1_3 = r1*r1*r1;
            r2_3 = r2*r2*r2;
            
            r1_5 = r1*r1*r1_3;
            r2_5 = r2*r2*r2_3;
            
            ux = ((1-mu)*(x(k)+mu))/(r1_3) + mu*(x(k)+mu-1)/(r2_3);
            uy =(1-mu)*y(k)/(r1_3)   + mu *y(k)/(r2_3);
            
            x(k+1) = (1/(1+h2))*( h*(xd(k)-y(k)) + h2*(yd(k)+x(k)) + x(k)*(1 + h2) + y(k)*(h3 + h) - h2*ux - h3*uy);
            y(k+1) = h*(yd(k)+x(k)) - h*x(k+1) + y(k)*(h2+1) -h2*uy;
            xd(k+1) = h*((y(k+1)-y(k))/h + x(k)) - h*ux + xd(k)-y(k) + y(k+1)+ h*u(1);
            yd(k+1) = h*((x(k) - x(k+1))/h + y(k)) - h*uy + yd(k)+x(k) - x(k+1) + h*u(2);
            
            % calculate the jacobian for the costate EOMs
            Uxx =  (1-mu) * ( 1/(r1_3) - ( 3*((x(k) + mu)^2)/(r1_5) ) ) + mu * (1/(r2_3) - 3*(x(k) - 1 + mu)^2 /(r2_5));
            Uyy = (1-mu) * ( 1/(r1_3) -  3*y(k)^2/(r1_5) ) + mu * ( 1/(r2_3) -3*y(k)^2/(r2_5));
            
            Uxy = -3*(1-mu)*y(k)*(x(k)+mu)/(r1_5) - 3*mu*y(k)*(x(k)-1+mu)/(r2_5); Uyx = Uxy;
            
            jacobian = zeros(4,4);
            % first row of jacobian is [df1/dx df1/dy df1/dxd df1/dyd] where f1 is the
            % update equation for x(k+1)
            
            jacobian(1,1) = 1/(1+h2) * (1+2*h2-h2*Uxx -h2*Uxy);
            jacobian(1,2) = 1/(1+h2) * (h3-h2*Uyx - h3*Uyy );
            jacobian(1,3) = h/(1+h2);
            jacobian(1,4) = h2/(1+h2);
            % second row of jacobian is [df2/dx df2/dy df2/dxd df2/dyd] where f2 is the
            % update equation for y(k+1)
            jacobian(2,1) = h -h*jacobian(1,1) - h2*Uyx;
            jacobian(2,2) = -h*jacobian(1,2) + h2+1 - h2*Uyy;
            jacobian(2,3) = -h*jacobian(1,3);
            jacobian(2,4) = h -h*jacobian(1,4);
            
            % third row of jacobian is [df3/dx df3/dy df3/dxd df3/dyd] where f3 is the
            % update equation for xd(k+1)
            jacobian(3,1) = 2*jacobian(2,1) +h - h*Uxx ;
            jacobian(3,2) = 2*jacobian(2,2) -2 -h*Uxy;
            jacobian(3,3) = 2*jacobian(2,3) + 1;
            jacobian(3,4) = 2*jacobian(2,4);
            
            % forth row of jacobian is [df4/dx df4/dy df4/dxd df4/dyd] where f4 is the
            % update equation for yd(k+1)
            jacobian(4,1) = 2-2*jacobian(1,1)-h*Uyx ;
            jacobian(4,2) = h-2*jacobian(1,2) - h*Uyy;
            jacobian(4,3) = -2*jacobian(2,3) ;
            jacobian(4,4) = 1-2*jacobian(2,4);
            
            % solve for costates
            % transpose the jacobian for column multiplication
            J = jacobian';
            B = costate(k,:);
            
            f1x = J(1,1);
            f2x = J(1,2);
            f3x = J(1,3);
            f4x = J(1,4);
            
            f1y = J(2,1);
            f2y = J(2,2);
            f3y = J(2,3);
            f4y = J(2,4);
            
            f1xd = J(3,1);
            f2xd = J(3,2);
            f3xd = J(3,3);
            f4xd = J(3,4);
            
            f1yd = J(4,1);
            f2yd = J(4,2);
            f3yd = J(4,3);
            f4yd = J(4,4);
            
            a = -f1y/f1x; b = -f1xd/f1x; c = -f1yd/f1x;
            e = -(f2xd+b*f2x)/(f2y+a*f2x);
            f = -(f2yd+c*f2x)/(f2y+a*f2x);
            g = - (f3yd + c*f3x + f*(f3y + a*f3x))/(f3xd + b*f3x + e*(f3y + a*f3x));
            
            a11 = f1x; a12 = f2x; a13 = f3x; a14 = f4x;
            a22 = f2y + a*f2x; a23 = f3y + a*f3x; a24 = f4y + a*f4x;
            a33 = f3xd + b*f3x + e*(f3y+a*f3x); a34 = f4xd + b*f4x + e*(f4y + a*f4x);
            a44 = f4yd + c*f4x + f*(f4y + a*f4x) + g*(f4xd + b*f4x + e*(f4y + a*f4x));
            
            ref = [a11 a12 a13 a14;...
                0 a22 a23 a24;...
                0 0 a33 a34; ...
                0 0 0 a44];
            
            beta1 = B(1);
            beta2 = B(2) + a*B(1);
            beta3 = B(3) + b*B(1) + e*(B(2) + a*B(1));
            beta4 = B(4) + c*B(1) + f*(B(2) + a*B(1)) + g*(B(3) + b* B(1) + e*(B(2) + a*B(1)));
            
            X = zeros(4,1);
            X(4) = beta4/a44;
            X(3) = beta3/a33 - a34/a33*X(4);
            X(2) = beta2/a22 - a23/a22*X(3) - a24/a22*X(4);
            X(1) = beta1/a11 - a12/a11*X(2) - a13/a11*X(3) - a14/a11*X(4);
            
            costate(k+1,:) = X';
            
            %             costate(k+1,:) = costate(k,:)/(jacobian);
            H(k) = costate(k+1,:)*state(k+1,:)';
        end
        
        state = [x y xd yd];
        
        
    case 'trap'
        
        x = zeros(num_steps,1);
        y = zeros(num_steps,1);
        xd = zeros(num_steps,1);
        yd = zeros(num_steps,1);
        
        Pxh = zeros(4,4,num_steps);
        Phh = zeros(4,4,num_steps);
        Pxx = zeros(4,4,num_steps);
        Phx = zeros(4,4,num_steps);
        
        x(1) = x0(1);
        y(1) = x0(2);
        xd(1) = x0(3);
        yd(1) = x0(4);
        
        Pxh(:,:,1) = Pxh0;
        Phh(:,:,1) = Phh0;
        Pxx(:,:,1) = Pxx0;
        Phx(:,:,1) = Phx0;
        
        for k = 1:num_steps-1
            switch constants.control_switch
                case 'off'
                    u = zeros(2,1);
                case 'on'
                    % express control vector in terms of costates
                    hv = costate(k,3:4);
                    u = -constants.um *hv/norm(hv);
                case 'sub'
                    rand_vec = constants.sub_opt;
                    u = constants.um*rand_vec/norm(rand_vec);
            end
            
            h2 = h*h;
            h3 = h*h*h;
            
            r1k = sqrt((x(k)+mu)^2 + y(k)^2 );
            r2k = sqrt( (x(k)+mu-1)^2+ y(k)^2 );
            
            r1k_3 = r1k*r1k*r1k;
            r2k_3 = r2k*r2k*r2k;
            
            r1k_5 = r1k_3*r1k*r1k;
            r2k_5 = r2k_3*r2k*r2k;
            
            uxk = ((1-mu)*(x(k)+mu))/(r1k_3) + mu*(x(k)+mu-1)/(r2k_3);
            uyk =(1-mu)*y(k)/r1k_3  + mu *y(k)/r2k_3;
            
            x(k+1) = 1/(1+h2) *(h*xd(k) -h*y(k) +h2*yd(k) + h2*x(k) + x(k)*(1+h2/2) + y(k)*(h+1/2*h3) ...
                -1/2*h3*uyk -1/2*h2*uxk);
            y(k+1) = h*yd(k) + h*x(k)-h*x(k+1) + y(k) + 1/2*h2*y(k)-h2/2*uyk;
            
            r1kp = sqrt((x(k+1)+mu)^2 + (y(k+1))^2);
            r2kp = sqrt((x(k+1)-1+mu)^2 + (y(k+1))^2);
            
            r1kp_3 = r1kp*r1kp*r1kp;
            r2kp_3 = r2kp*r2kp*r2kp;
            
            r1kp_5 = r1kp_3*r1kp*r1kp;
            r2kp_5 = r2kp_3*r2kp*r2kp;
            
            uxkp = ((1-mu)*(x(k+1)+mu))/(r1kp_3) + mu*(x(k+1)+mu-1)/(r2kp_3);
            uykp =(1-mu)*y(k+1)/r1kp_3  + mu *y(k+1)/r2kp_3;
            
            xd(k+1) = xd(k) - 2*y(k) + 2*y(k+1) + h/2*(x(k+1) + x(k)) - h/2*uxkp - h/2*uxk + h*u(1);
            yd(k+1) = yd(k) + 2*x(k) - 2*x(k+1) + h/2*(y(k+1) + y(k)) - h/2*uykp - h/2*uyk + h*u(2);
            
            % calculate the jacobian
            % second order partials at step k
            Uxxk =  (1-mu) * ( 1/r1k_3 - ( 3*((x(k) + mu)^2)/r1k_5 ) ) + mu * (1/r2k_3 - 3*(x(k) - 1 + mu)^2 /r2k_5);
            Uyyk = (1-mu) * ( 1/r1k_3 -  3*y(k)^2/r1k_5 ) + mu * ( 1/r2k_3 -3*y(k)^2/r2k_5);
            Uxyk = -3*(1-mu)*y(k)*(x(k)+mu)/r1k_5 - 3*mu*y(k)*(x(k)-1+mu)/r2k_5; Uyxk = Uxyk;
            
            jacobian = zeros(4,4);
            
            jacobian(1,1) = 1/(1+h2) *( h2 + 1 + 1/2*h2 -1/2*h3 * Uyxk - 1/2*h2 * Uxxk);
            jacobian(1,2) = 1/(1+h2) *(-h + h +h3/2 - h3/2*Uyyk - h2/2 * Uxyk );
            jacobian(1,3) = h/(1+h2);
            jacobian(1,4) = h2/(1+h2);
            
            jacobian(2,1) = h - h*jacobian(1,1) - h2/2 * Uyxk;
            jacobian(2,2) = -h*jacobian(1,2) + 1 + h2/2 -h2/2*Uyyk;
            jacobian(2,3) = -h*jacobian(1,3);
            jacobian(2,4) = h - h*jacobian(1,4);
            
            % partials of potential at k+1
            
            r1kp_partial = zeros(1,4);
            r1kp_partial(1) =  ((x(k+1) + mu) * jacobian(1,1) + y(k+1)*jacobian(2,1));
            r1kp_partial(2) =  ((x(k+1) + mu) * jacobian(1,2) + y(k+1)*jacobian(2,2));
            r1kp_partial(3) =  ((x(k+1) + mu) * jacobian(1,3) + y(k+1)*jacobian(2,3));
            r1kp_partial(4) =  ((x(k+1) + mu) * jacobian(1,4) + y(k+1)*jacobian(2,4));
            
            r2kp_partial = zeros(1,4);
            r2kp_partial(1) = ((x(k+1) - 1 + mu) * jacobian(1,1) + y(k+1)*jacobian(2,1));
            r2kp_partial(2) = ((x(k+1) - 1 + mu) * jacobian(1,2) + y(k+1)*jacobian(2,2));
            r2kp_partial(3) = ((x(k+1) - 1 + mu) * jacobian(1,3) + y(k+1)*jacobian(2,3));
            r2kp_partial(4) = ((x(k+1) - 1 + mu) * jacobian(1,4) + y(k+1)*jacobian(2,4));
            
            Uxxkp =  (1-mu) * (jacobian(1,1)/r1kp_3  - ( 3*((x(k+1) + mu))/r1kp_5 ) * r1kp_partial(1) ) ...
                + mu * (jacobian(1,1)/r2kp_3 - (3*(x(k+1) - 1 + mu)/r2kp_5) * r2kp_partial(1));
            Uxykp =  (1-mu) * (jacobian(1,2)/r1kp_3  - ( 3*((x(k+1) + mu))/r1kp_5 ) * r1kp_partial(2) ) ...
                + mu * (jacobian(1,2)/r2kp_3 - 3*(x(k+1) - 1 + mu)/r2kp_5 * r2kp_partial(2));
            Uxxdkp = (1-mu) * (jacobian(1,3)/r1kp_3 - ( 3*((x(k+1) + mu))/r1kp_5 ) * r1kp_partial(3) ) ...
                + mu * (jacobian(1,3)/r2kp_3 - (3*(x(k+1) - 1 + mu)/r2kp_5) * r2kp_partial(3));
            Uxydkp = (1-mu) * (jacobian(1,4)/r1kp_3 - ( 3*((x(k+1) + mu))/r1kp_5 ) * r1kp_partial(4) ) ...
                + mu * (jacobian(1,4)/r2kp_3 - 3*(x(k+1) - 1 + mu)/r2kp_5 * r2kp_partial(4));
            
            Uyxkp = (1-mu) * (jacobian(2,1)/r1kp_3 -  3*y(k+1)/r1kp_5 *r1kp_partial(1) ) ...
                + mu * ( jacobian(2,1)/r2kp_3 -3*y(k+1)/r2kp_5 * r2kp_partial(2));
            Uyykp = (1-mu) * (jacobian(2,2)/r1kp_3 -  3*y(k+1)/r1kp_5 *r1kp_partial(2) ) ...
                + mu * (jacobian(2,2)/r2kp_3 -3*y(k+1)/r2kp_5 * r2kp_partial(2));
            Uyxdkp = (1-mu) * (jacobian(2,3)/r1kp_3 -  3*y(k+1)/r1kp_5 *r1kp_partial(3) ) ...
                + mu * (jacobian(2,3)/r2kp_3 -3*y(k+1)/r2kp_5 * r2kp_partial(3));
            Uyydkp = (1-mu) * (jacobian(2,4)/r1kp_3 -  3*y(k+1)/r1kp_5 *r1kp_partial(4) ) ...
                + mu * (jacobian(2,4)/r2kp_3 -3*y(k+1)/r2kp_5 * r2kp_partial(4));
            
            % jacobian of velocity
            jacobian(3,1) = 2*jacobian(2,1) + h/2*(jacobian(1,1) + 1) - h/2*Uxxkp - h/2 * Uxxk;
            jacobian(3,2) = -2+ 2* jacobian(2,2) + h/2*jacobian(1,2) - h/2 * Uxykp - h/2 * Uxyk;
            jacobian(3,3) = 1+2*jacobian(2,3) + h/2*jacobian(1,3) - h/2*Uxxdkp;
            jacobian(3,4) = 2*jacobian(2,4) + h/2*jacobian(1,4) - h/2*Uxydkp;
            
            jacobian(4,1) = 2 - 2 * jacobian(1,1) + h/2*jacobian(2,1) -h/2*Uyxkp - h/2 * Uyxk;
            jacobian(4,2) = -2*jacobian(1,2) + h/2*(jacobian(2,2)+1) - h/2*Uyykp - h/2*Uyyk;
            jacobian(4,3) = -2*jacobian(1,3) + h/2 *jacobian(2,3) - h/2*Uyxdkp;
            jacobian(4,4) = 1 -2*jacobian(1,4) + h/2*jacobian(2,4) - h/2 * Uyydkp;
            
            % solve for costates
            % transpose the jacobian for column multiplication
            J = jacobian';
            B = costate(k,:);
            
            f1x = J(1,1);
            f2x = J(1,2);
            f3x = J(1,3);
            f4x = J(1,4);
            
            f1y = J(2,1);
            f2y = J(2,2);
            f3y = J(2,3);
            f4y = J(2,4);
            
            f1xd = J(3,1);
            f2xd = J(3,2);
            f3xd = J(3,3);
            f4xd = J(3,4);
            
            f1yd = J(4,1);
            f2yd = J(4,2);
            f3yd = J(4,3);
            f4yd = J(4,4);
            
            a = -f1y/f1x; b = -f1xd/f1x; c = -f1yd/f1x;
            e = -(f2xd+b*f2x)/(f2y+a*f2x);
            f = -(f2yd+c*f2x)/(f2y+a*f2x);
            g = - (f3yd + c*f3x + f*(f3y + a*f3x))/(f3xd + b*f3x + e*(f3y + a*f3x));
            
            a11 = f1x; a12 = f2x; a13 = f3x; a14 = f4x;
            a22 = f2y + a*f2x; a23 = f3y + a*f3x; a24 = f4y + a*f4x;
            a33 = f3xd + b*f3x + e*(f3y+a*f3x); a34 = f4xd + b*f4x + e*(f4y + a*f4x);
            a44 = f4yd + c*f4x + f*(f4y + a*f4x) + g*(f4xd + b*f4x + e*(f4y + a*f4x));
            
            ref = [a11 a12 a13 a14;...
                0 a22 a23 a24;...
                0 0 a33 a34; ...
                0 0 0 a44];
            
            beta1 = B(1);
            beta2 = B(2) + a*B(1);
            beta3 = B(3) + b*B(1) + e*(B(2) + a*B(1));
            beta4 = B(4) + c*B(1) + f*(B(2) + a*B(1)) + g*(B(3) + b* B(1) + e*(B(2) + a*B(1)));
            
            X = zeros(4,1);
            X(4) = beta4/a44;
            X(3) = beta3/a33 - a34/a33*X(4);
            X(2) = beta2/a22 - a23/a22*X(3) - a24/a22*X(4);
            X(1) = beta1/a11 - a12/a11*X(2) - a13/a11*X(3) - a14/a11*X(4);
            
            costate(k+1,:) = X';
            
%          costate(k+1,:) = costate(k,:)/(jacobian);
            H(k) = costate(k+1,:)*state(k+1,:)';
            
            % calculate the sensitivities
            % linearize about current state
%             cur_state = [x(k) y(k) xd(k) yd(k)];
%             cur_costate = costate(k,:);
%             [dfdx, dfdh, dgdx, dgdh, ~,~] = pcrtbp_trap_lin(cur_state,cur_costate,h,constants); % linearize system about current trajectory
%             
%             Pxh(:,:,k+1) = dfdx*Pxh(:,:,k)+dfdh*Phh(:,:,k);
%             Phh(:,:,k+1) = dgdx*Pxh(:,:,k)+dgdh*Phh(:,:,k);
%             Pxx(:,:,k+1) = dfdx*Pxx(:,:,k)+dfdh*Phx(:,:,k);
%             Phx(:,:,k+1) = dgdx*Pxx(:,:,k)+dgdh*Phx(:,:,k);
        end
        
        state = [x y xd yd];
        
end % end of mode switch

state = [state costate];
sens.Pxh = Pxh;
sens.Phh = Phh;
sens.Pxx = Pxx;
sens.Phx = Phx;
end % end of variational integrator function