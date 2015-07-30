function [Es,Eu,Ec,Vs,Vu,Vc]=pcrtbp_eig(eq_state, mu)

%   Author:
%       - Shankar Kulumani 4 September 2014
%           - corrected error when only inputting a PCRTBP state (4x1) no
%           need to modify the output of the Df_mat

% seperate out the eigenvector/values into the stable/unstable/center
% subspaces
tol = 1e-13;

[Df_mat] = linearized_eom_mat(eq_state, mu); % only need the upper 4x4 for the planar case


[V D ] = eig(Df_mat);  % V eigenvector - D eigenvalues

% remove any values that are too close to zero
small_real = abs(real(V)) < tol;
V(small_real) = imag(V(small_real))*1i;

small_imag = abs(imag(V)) < tol;
V(small_imag) = real(V(small_imag)); 

small_real = abs(real(D)) < tol;
D(small_real) = imag(D(small_real))*1i;

small_imag = abs(imag(D)) < tol;
D(small_imag) = real(D(small_imag)); 

% convert the square eigenvalue matrix into a 1D array
D = D(logical(eye(size(D))))';

% move all stable (left hand plane) eigenvalues/vectors
s_index = real(D) < 0 ; % stable
u_index = real(D) > 0 ; % unstable
c_index = abs(real(D)) < tol; % center/marginal

Es = D(s_index);
Eu = D(u_index);
Ec = D(c_index);

Vs = V(:,s_index);
Vu = V(:,u_index);
Vc = V(:,c_index);

