% 24 October 2014 - loads constants structure for various problems

function constants = crtbp_constants

% global ode parameters
constants.RelTol = 1e-9;
constants.AbsTol = 1e-9;
constants.ode_options = odeset('RelTol',constants.RelTol,'AbsTol',constants.AbsTol);
% global iteration tolerance
constants.tol = 1e-8; 

mu = 0.0125;
constants.l_scale = 384400; % distance of M2 from M1 (SMA) km
constants.t_scale = 2.361e6; % period of M2 about M1 seconds
% mu = 3.04036e-6;
% constants.l_scale = 1.496e8; % distance of M2 from M1 (SMA) km
% constants.t_scale = 3.147e7; % period of M2 about M1 seconds
constants.n_scale = 2*pi/constants.t_scale;% mean motion of M2 about M1 rad/sec
% vdim = v_scale * v_nondim
constants.v_scale = constants.n_scale*constants.l_scale; % mean velocity of M2 about M1 km/sec
% nondimensional to dimensional - adim = a_scale*a_nondim
constants.a_scale = constants.n_scale*constants.n_scale*constants.l_scale;
constants.sc_mass = 500; % kilogram mass of spacecraft
constants.km2meter = 1000/1; % convert kilometers to meters

constants.plot_scale = 'nondim'; % nondim or dim 
constants.mu = mu;
l_point = 1;
constants.l_point = l_point;
constants.Az = 500/constants.l_scale;
constants.earth_rad = 6378.137/constants.l_scale; % non-dimensionalized
% differential correction parameters
constants.diffcorr_select_pcrtbp = 'initial_ydot'; % initial_ydot or initial_x
constants.diffcorr_select_crtbp  = 'hold_z'; % hold_x or hold_z only for crtbp halo
constants.diffcorr_plot = 0;
constants.periodic_orbit_switch = 'halo'; % planar or halo for CRTBP case
constants.halo_switch = 'north';

% manifold parameters
constants.numsteps = 2000;
constants.manifold_steps = 200;
constants.manifold_tf_multipler = 5; % 
constants.epsilon = 50/constants.l_scale;
constants.manifold_plot = 0;

% poincare parameters
constants.poincare_steps = 1000;
constants.poincare_section = 2;

constants.filename = 'U3_L2_earthmoon.mat';

% pick energy level in case 3 (between energies of the L2 and L3 points to
% open up a neck in the secondary primary region
% calculate the energy of the L2 and L3 points
L_points = libration_points(mu);
L1 = [L_points(1,:) 0 0 0 0];
L2 = [L_points(2,:) 0 0 0 0];
L3 = [L_points(3,:) 0 0 0 0];

E1 = energyconst(L1,mu);
E2 = energyconst(L2,mu);
E3 = energyconst(L3,mu);

e_desired = E2 + 0.39*abs(E2-E3); % desired energy level of periodic orbit larger than E2 or more positive
% e_desired = -1.4778;

constants.e_desired = e_desired;
