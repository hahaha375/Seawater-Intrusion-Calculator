% This code used the Finited-Difference Methods (FDM) to solve the L
% equation. The description can be seen in the word file
% 3D unconfined aquifer head
clear;close all
%% ============================Here is the data set========================
% Define the grid size and the boundary functions
Lx = 100; % Length along x direction
Ly = 100; % Length along y direction
Lz = 10; % Length along z direction
m = 125+1; % number of grid points along x direction
n = 125+1; % number of grid points along y direction
N = 20+1; % number of grid points along z direction
delta_x = Lx/(m-1);
delta_y = Ly/(n-1);
H = 10; % be careful with the height
Hs = 9.4;
HL = 9.8;

% recharge term
I = 0.365/365*ones(n, m);

% well information
xw = 50; xw_nod = round(xw/delta_x)+1;
yw = 50; yw_nod = round(yw/delta_y)+1;
Qw = 25; % m3/d
I(yw_nod, xw_nod) = -Qw/delta_x/delta_y; % treat the pumping as recharge

% dispersion
alphaT = 0.4*0.1;
alpha = 1+0.025*(1-(alphaT/H)^(1/4));

% Processing stochastic field
x = linspace(0, Lx, m); % Grid coordinates in the x direction
y = linspace(0, Ly, n); % Grid coordinates in the y direction
z = linspace(0, Lz, N); % Grid coordinates in the z direction
[X1, Y1, Z1] = meshgrid(x, y, z);
K0 = load('K1.mat').K0; % Attention!!!! y,x,z
K0 = K0(:,:,end:-1:1);% Be careful with this command for unconfined
%aquifer!!!! The inversed matrix could increase the thickness of the bottom
%layer!!!! Be careful!!!
[X0, Y0, Z0] = meshgrid(linspace(0,Lx,size(K0,2)), linspace(0,Ly,size(K0,1)), linspace(0,Lz,size(K0,3)));
K0 = interp3(X0, Y0, Z0, K0, X1, Y1, Z1);

% iteration weighting factor
% recommendation: no vertical heterogeneity 1;
%                 vertical heterogeneity, 0.5~1
ssor = 0.8;
max_iter = 100;

%% ==========================End of the data set===========================
% Find the vertical mean
K = mean(K0, 3);
K = [K, K(:,end)];
K = [K; K(end,:)];
% Define the IC
Phi_0 = 0.5*alpha*Hs^2;
Phi_L = 0.5*HL^2;
phi1 = Hs*ones(n, m); % initial free surface
zeta = 0;
xt = 0.5*alpha^2*Hs^2; % potential at xt

%% ===============================caluculation=============================
% Each line represent a equation, Each colum represent a variable
% The reshape function prioritizes sorting by column, so we should numbered
% in column order
% ghost point method was used!!!!!!!
for iter = 1:max_iter
    fprintf('iteration %d \n', iter)
    % Attention!!!!!!!!!!!!!!!
    % Since the second loop, the dimensions of the matrix have changed,
    % so the coefficient matrix needs to be reset
    % Construct the coefficient matrix A and the right-hand side vector F
    index_i = zeros(m*n, 5);
    index_j = zeros(m*n, 5);
    value = zeros(m*n, 5);
    % A = zeros(m*n, m*n); % initialize A as a zero matrix
    F = zeros(m*n, 1); % initialize F as a zero vector
    zeta0 = zeta;
    phi0 = phi1;

    for i = 1:m % i represent the i index in the difference equation along x direction
        for j =1:n % j represent the j index in the difference equation along y direction
            k = (i-1)*n+j; % index of the (i, j) point
            %A(k, k) = -2/delta_x^2-2/delta_y^2; % set the diagonal element of A
            index_i(k, 1) = k;
            index_j(k, 1) = k;
            % Phi(i,j)
            value(k, 1) = -(K(j,i+1)-K(j,i))/delta_x^2-2*K(j,i)/delta_x^2-(K(j+1,i)-K(j,i))/delta_y^2-2*K(j,i)/delta_y^2;
            F(k) = -I(j, i);  % recharge term
            if i >1 % if not on the left boundary
                %A(k, k-n) = 1/delta_x^2; % set the left neighbor of A(k, k)
                index_i(k, 2) = k;
                index_j(k, 2) = k-n;
                % Phi(i-1,j)
                value(k, 2) = K(j,i)/delta_x^2;
            else % on the left boundary (constant head)
                F(k) = -K(j,i)*Phi_0/delta_x^2-I(j, i);
            end
            if i < m % if not on the right boundary
                %A(k, k+n) = 1/delta_x^2; % set the right neighbor of A(k, k)
                index_i(k, 3) = k;
                index_j(k, 3) = k+n;
                % Phi(i+1,j)
                value(k, 3) = 1*K(j,i+1)/delta_x^2;
            else % on the right boudnary (constan head)
                %A(k, k) = -1/delta_x^2-2/delta_y^2;
                index_i(k, 1) = k;
                index_j(k, 1) = k;
                F(k) = -K(j, i+1)*Phi_L/delta_x^2-I(j, i);
            end
            if j>1 % if not on the top boundary
                %A(k, k-1) = 1/delta_y^2;
                index_i(k, 4) = k;
                index_j(k, 4) = k-1;
                % Phi(i,j-1)
                value(k, 4) = K(j,i)/delta_y^2;
            else
                %A(k, k) = -2/delta_x^2-1/delta_y^2;
                index_i(k, 1) = k;
                index_j(k, 1) = k;
                value(k, 1) = -(K(j,i+1)-K(j,i))/delta_x^2-2*K(j,i)/delta_x^2-(K(j+1,i)-K(j,i))/delta_y^2-K(j,i)/delta_y^2;
            end
            if j<n %if not on the bottom boudnary
                %A(k, k+1) = 1/delta_y^2;
                index_i(k, 5) = k;
                index_j(k, 5) = k+1;
                % Phi(i,j+1)
                value(k, 5) = K(j+1,i)/delta_y^2;
            else
                %A(k, k) = -2/delta_x^2-1/delta_y^2;
                index_i(k, 1) = k;
                index_j(k, 1) = k;
                value(k, 1) = -(K(j,i+1)-K(j,i))/delta_x^2-2*K(j,i)/delta_x^2-K(j,i)/delta_y^2;
            end
        end
    end
    % for the four special point
    % only the point with two flux boundary is special
    % % right up
    % index_i((m-1)*n+1, 1) = (m-1)*n+1;
    % index_j((m-1)*n+1, 1) = (m-1)*n+1;
    % value((m-1)*n+1, 1) = -K(1,m)/delta_x^2-(K(2,m)-K(1,m))/delta_y^2-K(1,m)/delta_y^2;
    % % right down
    % index_i(m*n, 1) = m*n;
    % index_j(m*n, 1) = m*n;
    % value(m*n, 1) = -K(n,m)/delta_x^2-K(n,m)/delta_y^2;
    % Convert to a sparse matrix
    index_i = nonzeros(reshape(index_i, [], 1));
    index_j = nonzeros(reshape(index_j, [], 1));
    value = nonzeros(reshape(value, [], 1));
    A = sparse(index_i, index_j, value, m*n, m*n);
    F = sparse(F);
    % Solve the equation
    [L,U] = ilu (A); % use preconditioner, very useful
    tol = 1e-12;
    maxit = size(A,1);
    % select solver
    %Phi_star = gmres(A, F, [], tol, maxit, L, U, Phi_0*ones(maxit, 1));
    %Phi_star = bicgstab(A, F, tol, maxit, L, U, Phi_0*ones(maxit, 1));
    Phi_star = A\F;
    % change the shape
    Phi_star = reshape(Phi_star, n, m);
    % Solving head
    phi1 = zeros(size(Phi_star));
    Phi_star = full(Phi_star); % exclude the pumping well
    M_toe = contourc(x, y, Phi_star, [xt xt]);
    index = M_toe(2, 1); % seperate the toe locations
    x_toe = M_toe(1, 2:index);
    y_toe = M_toe(2, 2:index);
    [y_toe, sort_index] = sort(y_toe); % sort the toe location
    x_toe = x_toe(sort_index);
    y_toe = [0, y_toe, Ly, Ly, 0, 0]; % Create a closed polygon
    x_toe = [x_toe(1), x_toe, x_toe(end), 0, 0, x_toe(1)]; % use the polygon filtering the toe locations
    filter1 = inpolygon(squeeze(X1(:,:,1)),squeeze(Y1(:,:,1)),x_toe,y_toe) &...
                Phi_star<=xt; % (the potential will below the sea, leading to complex number)
    temp = Phi_star(filter1);
    phi1(filter1) = ((2*temp-alpha*Hs^2)*(1-1/alpha)).^0.5+Hs;
    phi1(filter1) = real(phi1(filter1))-imag(phi1(filter1));
    temp = Phi_star(~filter1);
    phi1(~filter1) = (2*temp).^0.5;
    % the freshwater surface cannot exceed the aquifer top
    filter2 = phi1>H;
    temp = Phi_star(filter2);
    phi1(filter2) = (temp+0.5*H^2)/H;

    % Solving interface
    zeta = alpha/(alpha-1)*Hs-1/(alpha-1)*phi1;
    zeta_temp = zeta;
    zeta(~filter1) = 0; % process the pumping well
    zeta(zeta>phi1) = phi1(zeta>phi1); % (the interface cannot exceed the free surface)
    nr1 = floor(zeta);
    nr2 = ceil(zeta);
    % calculate error
    Error = mean(zeta-zeta0, "all");
    Error2 = mean(phi1-phi0, "all");
    % jump out the iteration
    if abs(Error) < 1e-5 && abs(Error2) < 1e-5
        fprintf('Reach specified error at ietration %d with error %e and error2 %e \n', iter, Error, Error2)
        break
    end
    % Takes the average of the assumed value and the current calculated value
    phi1 = (1-ssor)*phi0+ssor*phi1;
    zeta = (1-ssor)*zeta0+ssor*zeta;

    % calculate meanK for next iteration method3
    interp_num = N+iter;
    for nn = 1:size(nr1,1)
        for mm = 1:size(nr1,2)
            if phi1(nn, mm)<H % process the free water surface
                z_interp = linspace(zeta(nn, mm), phi1(nn, mm), interp_num);
            else
                z_interp = linspace(zeta(nn, mm), H, interp_num);
            end
            K_interp = interp1(z, squeeze(K0(nn,mm,:)), z_interp);
            K(nn, mm) = (K_interp(1)/2+sum(K_interp(2:end-1))+K_interp(end)/2)/(interp_num-1);
        end
    end

    % processing BC
    K(:, end) = K(:, end-1);
    %K(1, :) = K(2, :);
    K(end, :) = K(end-1, :);

end


%============================plot the result===============================
% solution
figure()
[M1, ~] = contour(squeeze(X1(:,:,1)), squeeze(Y1(:,:,1)), Phi_star, [0.5*alpha^2*Hs^2 0.5*alpha^2*Hs^2]);
[M2, ~] = contour(squeeze(X1(:,:,1)), squeeze(Y1(:,:,1)), zeta_temp, [0 0]);
hold on
xlim([0 Lx])
ylim([0 Ly])
figure()
yy = zeta(51, :); % the middle slice
plot(x(yy>=0), yy(yy>=0), Color='k', LineWidth=4, LineStyle=':')
hold on
freesurface = phi1(51, :);
plot(x, freesurface, Color='b', LineWidth=4, LineStyle=':')
figure()
surf(phi1, 'EdgeColor','none');light;

% save the result
save('ZETA.mat', "zeta", "M1");