% This code used the Finited-Difference Methods (FDM) to solve the L
% equation. The description can be seen in the word file
clear;close all
tic
% confined aquifer !!!!
% constant-flux
% Define the grid size and the boundary functions
Lx = 100; % Length along x direction
Lz = 10; % Length along z direction
n = 400+1; % number of grid points along x direction
N = 40+1; % number of grid points along z direction
delta_x = Lx/(n-1);
H = 10;
Hs = 12;
Q0 = -0.6;
% recharge term
I = 0.365/365*ones(n, 1);
% well information
xw = 50; xw_nod = round(xw/delta_x)+1;
Qw = 0; % m3/d
I(xw_nod) = -Qw/delta_x; % treat the pumping as recharge

alphaT = 0.*0.1;
alpha = 1+0.025*(1-(alphaT/H)^(1/4));
% this is the end of the input

Phi_0 = alpha*H*Hs-0.5*alpha*H^2;
x = linspace(0, Lx, n); % Grid coordinates in the x direction
z = linspace(0, Lz, N); % Grid coordinates in the z direction
[X, Z] = meshgrid(x, z);

% Read Gaussian random fields
K0 = load('K1.mat').K;
K0 = K0(end:-1:1,:);
x0 = linspace(0, Lx, 400+1);
z0 = linspace(0, Lz, 40+1);
[X0, Z0] = meshgrid(x0, z0);
K0 = interp2(X0, Z0, K0, X, Z);

% iteration weighting factor
% recommendation: no vertical heterogeneity 1;
%                 vertical heterogeneity, 0.5~1
ssor = 0.8;

% Find the vertical mean
K = mean(K0, 1);
K = [K,K(end)]; % Add an identical element at the end of K to facilitate interpolation
% Define the IC
phi0 = (alpha*Hs-(alpha-1)*H)*ones(n, 1);
zeta = 0;
%K = 6*ones(size(K)); % Warning: The maximum toe position cannot exceed the model length!!!!!!
% Construct the coefficient matrix A and the right-hand side vector F
A = zeros(n, n); % initialize A as a zero matrix
F = zeros(n, 1); % initialize F as a zero vector
% Each line represent a equation, Each colum represent a variable
% The reshape function prioritizes sorting by column, so we should numbered
% in column order
% ghost point method was used!

for iter =1:100
    fprintf('iteration %d \n', iter)
    zeta0 = zeta;
    for i = 1:n % i represent the i index in the difference equation along x direction
        k = i; % index of the (i, j) point
        A(k, k) = -(K(i+1)-K(i))/delta_x^2-K(i)*2/delta_x^2; % set the diagonal element of A
        F(k) = -I(i);  % recharge term
        if i >1 % if not on the left boundary
            A(k, k-1) = K(i)*1/delta_x^2; % set the left neighbor of A(k, k)
        else % on the left boundary (constant head)
            F(k) = -K(i)*Phi_0/delta_x^2-I(i);
        end
        if i < n % if not on the right boundary
            A(k, k+1) = (K(i+1)-K(i))/delta_x^2+K(i)/delta_x^2; % set the right neighbor of A(k, k)
        else % on the right boudnary (constan flux)
            A(k, k) = -K(i)/delta_x^2;
            F(k) = Q0/delta_x-I(i);
        end
    end
    A = sparse(A);
    F = sparse(F);
    % Solve the equation
    % [L,U] = ilu (A); % use preconditioner, very useful
    % tol = 1e-7;
    % maxit = size(A,1);
    % Phi_star = gmres(A, F, [], tol, [], L, U, Phi_0*ones(maxit, 1));
    Phi_star = A\F;
    % Solving head
    phi1 = zeros(size(Phi_star));
    xt = alpha*H*Hs-0.5*H^2;
    temp = Phi_star(Phi_star<=xt);
    phi1(Phi_star<=xt) = ((temp-alpha*H*Hs+0.5*alpha*H^2)*2*(alpha-1)).^0.5-(alpha-1)*H+alpha*Hs;
    temp = Phi_star(Phi_star>xt);
    phi1(Phi_star>xt) = (temp+0.5*H^2)/H;
    % Solving interface
    zeta = alpha/(alpha-1)*Hs-1/(alpha-1)*phi1;
    zeta(zeta<0) = 0;
    nr1 = floor(zeta);
    nr2 = ceil(zeta);
    % calculate error
    Error = mean(zeta-zeta0, "all");
    % jump out the iteration
    if abs(Error) < 1e-5
        fprintf('Reach specified error at ietration %d with error %e \n', iter, Error)
        break
    end
    % Takes the average of the assumed value and the current calculated value
    phi1 = (1-ssor)*phi0+ssor*phi1;
    zeta = (1-ssor)*zeta0+ssor*zeta;
    % calculate meanK for next iteration
    interp_num = N+iter;
    for nn = 1:size(nr1,1)
        z_interp = linspace(zeta(nn), H, interp_num);
        K_interp = interp1(z, squeeze(K0(:, nn)), z_interp);
        K(nn) = (K_interp(1)/2+sum(K_interp(2:end-1))+K_interp(end)/2)/(interp_num-1);
    end
    K(end) = K(end-1);

end

elapsedTime = toc;
fprintf('The program ran for %.4f seconds.\n', elapsedTime);

plot(zeta); % plot the interface
hold on
plot(phi1, '--b') % plot the free water-surface