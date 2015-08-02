% 25 Feb 2015 - Rectangle rule for PCRTBP variational integrator

function [t, state] = pcrtbp_variational(x0,t0,tf, mode, num_steps, constants)

% initial state is defined as x0 = [x0 y0 xd0 yd0]
% transform from pos,vel to pos, mom
t = linspace(t0,tf,num_steps)';
h = t(2)-t(1);
mu = constants.mu;

u = zeros(2,1);

% switch integration scheme based on mode (rect, mid, trap)
switch mode
    case 'rect'
        % intialize arrays for the propogated states
        x = zeros(num_steps,1);
        y = zeros(num_steps,1);
        xd = zeros(num_steps,1);
        yd = zeros(num_steps,1);
        
        x(1) = x0(1,1);
        y(1) = x0(1,2);
        xd(1) = x0(1,3);
        yd(1) = x0(1,4);
        
        %         [x0] = mom_map(x0); % map from x,y,xd,yd to x,y,px,py
        %         % intialize arrays for the propogated states
        %         x = zeros(num_steps,1);
        %         y = zeros(num_steps,1);
        %         px = zeros(num_steps,1);
        %         py = zeros(num_steps,1);
        %
        %         x(1) = x0(1,1);
        %         y(1) = x0(1,2);
        %         px(1) = x0(1,3);
        %         py(1) = x0(1,4);
        
        for k = 1:num_steps-1
            %             r1 = sqrt((x(k)+mu)^2 + y(k)^2 );
            %             r2 = sqrt( (x(k)+mu-1)^2+ y(k)^2 );
            %
            %             A = h*(1-mu)*(x(k)+mu)/r1^3;
            %             B = h*mu*(x(k)-1+mu)/r2^3;
            %             C = h*(1-mu)*y(k)/r1^3;
            %             D = h*mu*y(k)/r2^3;
            %
            %             x(k+1) = (1/(1+h^2))*( h*px(k) + h^2*py(k) + x(k) + h^2 * x(k) + h^3 * y(k) + h*y(k) - h*A-h*B-h^2*C - h^2*D  );
            %             y(k+1) = h*py(k) - h*x(k+1) + h^2*y(k) + y(k) -h*C -h*D;
            %             px(k+1) = h*((y(k+1)-y(k))/h + x(k)) - A -B + px(k);
            %             py(k+1) = h*((x(k) - x(k+1))/h + y(k)) - C -D + py(k);
            
            r1 = sqrt((x(k)+mu)^2 + y(k)^2 );
            r2 = sqrt( (x(k)+mu-1)^2+ y(k)^2 );
            
            ux = ((1-mu)*(x(k)+mu))/(r1^3) + mu*(x(k)+mu-1)/(r2^3);
            uy =(1-mu)*y(k)/r1^3  + mu *y(k)/r2^3;
            
            x(k+1) = (1/(1+h*h))*( h*(xd(k)-y(k)) + h*h*(yd(k)+x(k)) + x(k)*(1 + h*h) + y(k)*(h*h*h + h) - h*h*ux - h*h*h*uy);
            y(k+1) = h*(yd(k)+x(k)) - h*x(k+1) + y(k)*(h*h+1) -h*h*uy;
            xd(k+1) = h*((y(k+1)-y(k))/h + x(k)) - h*ux + xd(k)-y(k) + y(k+1)+ h*u(1);
            yd(k+1) = h*((x(k) - x(k+1))/h + y(k)) - h*uy + yd(k)+x(k) - x(k+1) + h*u(2);
        end
        % map from x,y,px,py to x,y,xd,yd (inverse mom map)
        %         [state] = inv_mom_map([x y px py]);
        state = [x y xd yd];
    case 'trap'
%         [x0] = mom_map(x0); % map from x,y,xd,yd to x,y,px,py
%         % intialize arrays for the propogated states
%         x = zeros(num_steps,1);
%         y = zeros(num_steps,1);
%         px = zeros(num_steps,1);
%         py = zeros(num_steps,1);
%         
%         x(1) = x0(1,1);
%         y(1) = x0(1,2);
%         px(1) = x0(1,3);
%         py(1) = x0(1,4);
        % trapezoidal mode
                % intialize arrays for the propogated states
        x = zeros(num_steps,1);
        y = zeros(num_steps,1);
        xd = zeros(num_steps,1);
        yd = zeros(num_steps,1);
        
        x(1) = x0(1,1);
        y(1) = x0(1,2);
        xd(1) = x0(1,3);
        yd(1) = x0(1,4);
        
        for k = 1:num_steps-1
%             r1k = sqrt((x(k)+mu)^2 + (y(k))^2);
%             r2k = sqrt((x(k)-1+mu)^2 + (y(k))^2);
%             
%             A = (1-mu)*(x(k)+mu)/r1k^3;
%             B = mu*(x(k)-1+mu)/r2k^3;
%             C = (1-mu)*y(k)/r1k^3;
%             D = mu*y(k)/r2k^3;
%             
%             x(k+1) = 1/(1+h^2) *(h*px(k) +h^2*py(k) + x(k)*(1+h^2/2) + y(k)*(h+1/2*h^3) ...
%                 -1/2*h^3*(C+D) -1/2*h^2*(A+B));
%             y(k+1) = h*py(k)-h*x(k+1) + y(k) + 1/2*h^2*y(k)-h^2/2*(C+D);
%             
%             r1kp = sqrt((x(k+1)+mu)^2 + (y(k+1))^2);
%             r2kp = sqrt((x(k+1)-1+mu)^2 + (y(k+1))^2);
%             
%             Ap = (1-mu)*(x(k+1)+mu)/r1kp^3;
%             Bp = mu*(x(k+1)-1+mu)/r2kp^3;
%             Cp = (1-mu)*y(k+1)/r1kp^3;
%             Dp = mu*y(k+1)/r2kp^3;
%             
%             
%             
%             %     both answers should be the same
%             px(k+1) = px(k) + h/2*(x(k+1)+x(k)) + y(k+1) - y(k) - h/2*(Ap+Bp+A+B);
%             py(k+1) = py(k) + h/2*(y(k+1)+y(k)) + x(k) - x(k+1) - h/2*(Cp + Dp +C+D);
            
            %     px(k+1) = 1/2*((x(k+1)-x(k))/h-y(k)) + 1/2*((x(k+1)-x(k))/h-y(k+1)) + h/2*((y(k+1)-y(k))/h+x(k+1)) ...
            %         -h/2*(Ap+Bp);
            %     py(k+1) = 1/2*((y(k+1)-y(k))/h+x(k)) - h/2*((x(k+1)-x(k))/h-y(k+1)) + 1/2*((y(k+1)-y(k))/h+x(k+1)) ...
            %         -h/2*(Cp+Dp);
            
            %         px(k+1) = x(k+1)/h - x(k)/h -y(k) + h*x(k+1)/2 -h/2*(Ap + Bp);
            %         py(k+1) = y(k+1)/h - y(k)/h +x(k) + h/2*y(k+1) - h/2*(Cp+Dp);
            
            
            r1k = sqrt((x(k)+mu)^2 + y(k)^2 );
            r2k = sqrt( (x(k)+mu-1)^2+ y(k)^2 );
            
            uxk = ((1-mu)*(x(k)+mu))/(r1k^3) + mu*(x(k)+mu-1)/(r2k^3);
            uyk =(1-mu)*y(k)/r1k^3  + mu *y(k)/r2k^3;
            
            x(k+1) = 1/(1+h^2) *(h*xd(k) -h*y(k) +h^2*yd(k) + h^2*x(k) + x(k)*(1+h^2/2) + y(k)*(h+1/2*h^3) ...
                -1/2*h^3*uyk -1/2*h^2*uxk);
            y(k+1) = h*yd(k) + h*x(k)-h*x(k+1) + y(k) + 1/2*h^2*y(k)-h^2/2*uyk;
            
            r1kp = sqrt((x(k+1)+mu)^2 + (y(k+1))^2);
            r2kp = sqrt((x(k+1)-1+mu)^2 + (y(k+1))^2);
            
            uxkp = ((1-mu)*(x(k+1)+mu))/(r1kp^3) + mu*(x(k+1)+mu-1)/(r2kp^3);
            uykp =(1-mu)*y(k+1)/r1kp^3  + mu *y(k+1)/r2kp^3;
            
            xd(k+1) = xd(k) - 2*y(k) + 2*y(k+1) + h/2*(x(k+1) + x(k)) - h/2*uxkp - h/2*uxk + h*u(1);
            yd(k+1) = yd(k) + 2*x(k) - 2*x(k+1) + h/2*(y(k+1) + y(k)) - h/2*uykp - h/2*uyk + h*u(2);
            
        end
        
        % map from x,y,px,py to x,y,xd,yd (inverse mom map)
%         [state] = inv_mom_map([x y px py]);
state = [x y xd yd];
    otherwise
        fprintf('Wrong integration mode. Use rect or trap \n')
        
end


