% 21 July 2015
% interpolate between points on the Poincare section

% extract out the target states that are to the left and below
left_cross = cross_state(cross_state(:,1)< (1-constants.mu),:);
ld_cross = left_cross(left_cross(:,3) < 0,:);
[y,index] = unique(ld_cross(:,1));
ld_cross = ld_cross(index,:);

target_x = ld_cross(:,1);
target_xd = ld_cross(:,3);
% find the limits in the min x and max x direction
min_x = target_x(1);
max_x = target_x(end);

grid_target = linspace(min_x,max_x,100);
% use pchip to interpolate the target orbits
p_target = pchip(target_x, target_xd, grid_target);

figure(3)
hold all
plot(target_x, target_xd,'o',grid_target, p_target,'-')
% interpolate the reachable set now
[y,index] = unique(reach_poincare(:,1));
initial_poincare = reach_poincare(index,:);

initial_x = initial_poincare(:,1);
initial_xd = initial_poincare(:,3);
min_x = initial_x(1);
max_x = initial_x(end);
grid_initial = linspace(min_x,max_x,100);
p_initial = pchip(initial_x, initial_xd, grid_initial);

% plot the target poincare section and the interpolated values

figure(3)
hold all
plot(initial_x, initial_xd,'o',grid_initial,p_initial,'-')

% [x_int,xd_int,iout,jout] = intersections(grid_target,p_target,grid_initial,p_initial);
[x_int,xd_int,iout,jout] = intersections(target_x,target_xd,initial_x,initial_xd);
plot(x_int,xd_int,'g*','markersize',10)

% calculate the target state to control towards
[ydot] = energyfcn(x_int, 0, xd_int, constants.mu, energyconst(moon_x0',constants.mu));