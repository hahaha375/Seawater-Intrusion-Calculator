clear; close all
zeta = load('zeta.mat').zeta;
phi1 = load('zeta.mat').phi1;

% grid for FDM_SWI
Lx = 100; % Length along x direction
Ly = 100; % Length along y direction
Lz = 10;
m = 125+1; % number of grid points along x direction
n = 125+1; % number of grid points along y direction
delta_x = Lx/(m-1);
delta_y = Ly/(n-1);
x_mat = linspace(0, Lx, m); % Grid coordinates in the x direction
y_mat = linspace(0, Ly, n); % Grid coordinates in the y direction
[X_mat, Y_mat] = meshgrid(x_mat, y_mat);

% grid for SEAWAT
nx = 125;
ny = 125;
nz = 10;
dx = Lx/nx;
dy = Ly/ny;
dz = Lz/nz;
x_sea = dx/2:dx:Lx-dx/2;
y_sea = -Ly/2+dy/2:dy:Ly/2-dy/2; 
z_sea = dz/2:dz:Lz-dz/2; % Here we check whether the sharp interface is above the grid center point; this initial value is still reasonable
[X_sea, Y_sea] = meshgrid(x_sea, y_sea);

% head data for SEAWAT
Head_SEAWAT = ones(nz, ny, nx);
Head_temp = interp2(X_mat, Y_mat-Ly/2, phi1, X_sea, Y_sea);
for i = 1:nx
    for j = 1:ny
        Head_SEAWAT(:, j, i) = Head_temp(ny+1-j, i);
    end
end

save('Head_SEAWAT.mat', 'Head_SEAWAT')

% concentration data for SEAWAT
Con_SEAWAT = zeros(nz, ny, nx);
ZETA_temp = interp2(X_mat, Y_mat-Ly/2, zeta, X_sea, Y_sea);
for i = 1:nx
    for j = 1:ny
        interface = ceil(ZETA_temp(ny+1-j, i));
        if interface > 0 
            index_z = sum(z_sea>interface);
            Con_SEAWAT(1:index_z, j, i) = 0; % the grid number for seawater is from top to bottom
            Con_SEAWAT(index_z+1:end, j, i) = 35;
        else
            Con_SEAWAT(:, j, i) = 0;
        end
    end
end

save('Con_SEAWAT.mat', 'Con_SEAWAT')