%% Definition of Parameters
clear, clc;

% Domain Related
N_x = 120; % num of x nodes
%6451
N_y = N_x; % num of y nodes
d_x = 1; % dist between x nodes
d_y = 1; % dist between y nodes

% LBM Related
% Lattice velocity
ksi = [0 1 0 -1 0 1 -1 -1 1; ...
       0 0 1 0 -1 1 1 -1 -1 ];
min_error = 0.000001; % error when sim ends
w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]; % weights for D2Q9

c_s = 1/sqrt(3); % speed of sound (D2Q9)

Tau = 0.54; % relaxation time
Rho_in = 2;
vis = (Tau-0.5) * c_s^2; % kinematic viscosity

L = N_x; % Length of box
min_error = 0.000001; % error when sim ends

%% Initialization
Rho_ref=2;
Rho = ones(1, N_y, N_x)*Rho_ref; % Density
U = zeros(2, N_y, N_x); % Velocity

f = zeros(9, N_y, N_x); % PDF for all 9 directions at all locations

f = eqm_d2q9(Rho, U, ksi, w);

m = zeros(9, N_y, N_x);
m_eq = m;

f_new = f; % Update variable
f_eq = f; % Equilibrium

timer = 0;
max_timer = 10000;
cont = true;
%figure
%r = animatedline;
%title("Residuals")

% M = [ 1,  1,  1,  1,  1,  1,  1,  1,  1; ... 
%      -4, -1, -1, -1, -1,  2,  2,  2,  2; ...
%       4, -2, -2, -2, -2,  1,  1,  1,  1; ... 
%       0,  1,  0, -1,  0,  1, -1, -1,  1; ... 
%       0, -2,  0,  2,  0,  1, -1, -1,  1; ... 
%       0,  0,  1,  0, -1,  1,  1, -1, -1; ... 
%       0,  0, -2,  0,  2,  1,  1, -1, -1; ... 
%       0,  1, -1,  1, -1,  0,  0,  0,  0; ... 
%       0,  0,  0,  0,  0,  1, -1,  1, -1; ... 
%      ]; % transformation matrix

M=[1   1    1   1   1   1   1   1   1;...
    -4  -1  -1  -1  -1   2   2   2   2;...
    4   -2  -2  -2  -2   1   1   1   1;...
    0    1   0  -1    0   1  -1  -1  1;...
    0  -2    0   2    0   1  -1  -1  1;...
    0    0   1   0   -1   1   1  -1  -1;...
    0    0  -2   0    2   1   1  -1  -1;...
    0    1  -1   1   -1   0   0   0   0;...
    0    0   0   0     0   1  -1   1 -1];
%Mew = 0.015; % Kinematic viscosity
%U_lid = Re*Mew/L; % Velocity of the moving lid
%Tau = 3*Mew+0.5; % Relaxation Time

%U_top = Re*Mew/N_x; % Top velocity

Re = 100; % Reanolds number
Tau = 0.8;
vis = (Tau-0.5)*c_s^2;
U_top = Re*vis/N_x;
s_89 = 1/Tau;
s_2 = 0.001;
s_3 = 0.001;
s_57 = 0.001;
S = diag([0; s_2; s_3; 0; s_57; 0; s_57; s_89; s_89]); % relaxation vector

%% Solving
tic
res_list = zeros(1, max_timer);
super_guh = zeros(9, N_y, N_x); % temp
% Interior nodes
int_x = 2:(N_x-1);
int_y = 2:(N_x-1);

for t = 1:max_timer
    % Streaming / Boundary Conditions
    % Top left
    f_new(1, 1, 1) = f(1, 1  , 1  );
    f_new(3, 1, 1) = f(3, 1+1, 1  );
    f_new(4, 1, 1) = f(4, 1  , 1+1);
    f_new(7, 1, 1) = f(7, 1+1, 1+1);
    % Abnormal
    Rho_c = Rho(1, 1, 1+1);
    f_new(2, 1, 1) = f_new(4, 1, 1);
    f_new(5, 1, 1) = f_new(3, 1, 1);
    f_new(9, 1, 1) = f_new(7, 1, 1);
    f_new(6, 1, 1) = (Rho_c - f_new(1, 1, 1) - f_new(2, 1, 1) - f_new(3, 1, 1) - f_new(4, 1, 1) - f_new(5, 1, 1) - f_new(7, 1, 1) - f_new(9, 1, 1))/2;
    f_new(8, 1, 1) = f_new(6, 1, 1);

    % Top right
    f_new(1, 1, N_x) = f(1, 1  , N_x  );
    f_new(2, 1, N_x) = f(2, 1  , N_x-1);
    f_new(3, 1, N_x) = f(3, 1+1, N_x  );
    f_new(6, 1, N_x) = f(6, 1+1, N_x-1);
    % Abnormal
    Rho_c = Rho(1, 1, N_x-1);
    f_new(4, 1, N_x) = f_new(2, 1, N_x);
    f_new(5, 1, N_x) = f_new(3, 1, N_x);
    f_new(8, 1, N_x) = f_new(6, 1, N_x);
    f_new(1, 1, N_x) = (Rho_c - f_new(1, 1, N_x) - f_new(2, 1, N_x) - f_new(3, 1, N_x) - f_new(6, 1, N_x) - f_new(5, 1, N_x) - f_new(8, 1, N_x) - f_new(4, 1, N_x))/2;
    f_new(9, 1, N_x) = f_new(7, 1, N_x);

    % Top
    f_new(1, 1, int_x) = f(1, 1  , int_x  );
    f_new(2, 1, int_x) = f(2, 1  , int_x-1);
    f_new(3, 1, int_x) = f(3, 1+1, int_x  );
    f_new(4, 1, int_x) = f(4, 1  , int_x+1);
    f_new(6, 1, int_x) = f(6, 1+1, int_x-1);
    f_new(7, 1, int_x) = f(7, 1+1, int_x+1);
    % Abnormal
    p = f_new(1, 1, int_x) + f_new(2, 1, int_x) + f_new(4, 1, int_x) + 2*f_new(3, 1, int_x) + f_new(6, 1, int_x) + f_new(7, 1, int_x);
    f_new(5, 1, int_x) = f_new(3, 1, int_x); % Bounceback
    f_new(8, 1, int_x) = (-p*U_top + f_new(2, 1, int_x) + 2*f_new(6, 1, int_x) - f_new(4, 1, int_x))/2;
    f_new(9, 1, int_x) = ( p*U_top - f_new(2, 1, int_x) + 2*f_new(7, 1, int_x) + f_new(4, 1, int_x))/2;

    % Bottom left corner
    f_new(1, N_y, 1) = f(1, N_y  , 1  );
    f_new(4, N_y, 1) = f(4, N_y  , 1+1);
    f_new(5, N_y, 1) = f(5, N_y-1, 1  );
    f_new(8, N_y, 1) = f(8, N_y-1, 1+1);
    % Abnormal
    f_new(2, N_y, 1) = f(4, N_y, 1); % Bounceback
    f_new(3, N_y, 1) = f(5, N_y, 1); % Bounceback
    f_new(6, N_y, 1) = f(8, N_y, 1); % Bounceback
    f_new(7, N_y, 1) = (Rho_c - f_new(1, N_y, 1) - f_new(4, N_y, 1) - f_new(5, N_y, 1) - f_new(8, N_y, 1) - f_new(2, N_y, 1) - f_new(3, N_y, 1) - f_new(6, N_y, 1))/2;
    f_new(9, N_y, 1) = f_new(7, N_y, 1);

    % Bottom right corner
    f_new(1, N_y, N_x) = f(1, N_y  , N_x  );
    f_new(2, N_y, N_x) = f(2, N_y  , N_x-1);
    f_new(5, N_y, N_x) = f(5, N_y-1, N_x  );
    f_new(9, N_y, N_x) = f(9, N_y-1, N_x-1);
    % Abnormal
    Rho_c = Rho(1, N_x, N_x-1); % Approximation for rho
    f_new(3, N_y, N_x) = f_new(5, N_y, N_x); % Bounceback
    f_new(4, N_y, N_x) = f_new(2, N_y, N_x); % Bounceback
    f_new(7, N_y, N_x) = f_new(9, N_y, N_x); % Bounceback
    f_new(6, N_y, N_x) = (Rho_c - f_new(1, N_y, N_x) - f_new(2, N_y, N_x) - f_new(5, N_y, N_x) - f_new(9, N_y, N_x) - f_new(3, N_y, N_x) - f_new(4, N_y, N_x) - f_new(7, N_y, N_x))/2;
    f_new(8, N_y, N_x) = f_new(6, N_y, N_x);

    % Bottom boundary
    f_new(1, N_y, int_x) = f(1, N_y  , int_x  );
    f_new(2, N_y, int_x) = f(2, N_y  , int_x-1);
    f_new(4, N_y, int_x) = f(4, N_y  , int_x+1);
    f_new(5, N_y, int_x) = f(5, N_y-1, int_x  );
    f_new(8, N_y, int_x) = f(8, N_y-1, int_x+1);
    f_new(9, N_y, int_x) = f(9, N_y-1, int_x-1);
    % Abnormal
    f_new(3, N_y, int_x) = f_new(5, N_y, int_x); % Bounceback
    f_new(6, N_y, int_x) = f_new(8, N_y, int_x) + (f_new(4, N_y, int_x) - f_new(2, N_y, int_x))/2;
    f_new(7, N_y, int_x) = f_new(9, N_y, int_x) - (f_new(4, N_y, int_x) - f_new(2, N_y, int_x))/2;

    % Left boundary
    f_new(1, int_y, 1) = f(1, int_y  , 1  );
    f_new(3, int_y, 1) = f(3, int_y+1, 1  );
    f_new(4, int_y, 1) = f(4, int_y  , 1+1);
    f_new(5, int_y, 1) = f(5, int_y-1, 1  );
    f_new(7, int_y, 1) = f(7, int_y+1, 1+1);
    f_new(8, int_y, 1) = f(8, int_y-1, 1+1);
    % Abnormal
    f_new(2, int_y, 1) = f_new(4, int_y, 1); % Bounceback
    f_new(6, int_y, 1) = f_new(8, int_y, 1) + (f_new(5, int_y, 1) - f_new(3, int_y, 1))/2;
    f_new(9, int_y, 1) = f_new(7, int_y, 1) - (f_new(5, int_y, 1) - f_new(3, int_y, 1))/2;

    % Right boundary
    f_new(1, int_y, N_x) = f(1, int_y  , N_x  );
    f_new(2, int_y, N_x) = f(2, int_y  , N_x-1);
    f_new(3, int_y, N_x) = f(3, int_y+1, N_x  );
    f_new(5, int_y, N_x) = f(5, int_y-1, N_x  );
    f_new(6, int_y, N_x) = f(6, int_y+1, N_x-1);
    f_new(9, int_y, N_x) = f(9, int_y-1, N_x-1);
    % Abnormal
    f_new(4, int_y, N_x) = f_new(2, int_y, N_x); % Bounceback
    f_new(7, int_y, N_x) = f_new(9, int_y, N_x) + (f_new(5, int_y, N_x) - f_new(3, int_y, N_x))/2;
    f_new(8, int_y, N_x) = f_new(6, int_y, N_x) - (f_new(5, int_y, N_x) - f_new(3, int_y, N_x))/2;

    % All interior nodes
    f_new(1, int_y, int_x) = f(1, int_y  , int_x  );
    f_new(2, int_y, int_x) = f(2, int_y  , int_x-1);
    f_new(3, int_y, int_x) = f(3, int_y+1, int_x  );
    f_new(4, int_y, int_x) = f(4, int_y  , int_x+1);
    f_new(5, int_y, int_x) = f(5, int_y-1, int_x  );
    f_new(6, int_y, int_x) = f(6, int_y+1, int_x-1);
    f_new(7, int_y, int_x) = f(7, int_y+1, int_x+1);
    f_new(8, int_y, int_x) = f(8, int_y-1, int_x+1);
    f_new(9, int_y, int_x) = f(9, int_y-1, int_x-1);

    % % Collision
    % Rho_old = Rho;
    % % Rho, U calculation
    % [Rho, U] = rhoNu(f_new, ksi);
    % % f_eq calculation
    % f_eq = eqm_d2q9(Rho, U, ksi, w);
    % 
    % % BGK Collision and Update
    % f = f_new - (f_new-f_eq)/Tau;
    % 
    % [guh, res_list(t)] = res(Rho_old, Rho, min_error);
    
    %addpoints(r, t, max_error)
    %drawnow

    %progress(timer, t);
    Rho_old = Rho;
    Rho = sum(f_new, 1);
    U = pagemtimes(ksi,f_new)./Rho;

    % m = pagemtimes(M, f_new);
    % Rho = m(1, :, :);
    % U(1,:,:) = m(4,:,:);
    % U(2,:,:) = m(6,:,:);

    % super_guh(1,:,:) = ones(1, N_y, N_x);
    % super_guh(2,:,:) = -2 + 3*sum(U.^2, 1);
    % super_guh(3,:,:) = 1 - 3*sum(U.^2, 1);
    % super_guh(4,:,:) = U(1,:,:);
    % super_guh(5,:,:) = -U(1,:,:);
    % super_guh(6,:,:) = U(2,:,:);
    % super_guh(7,:,:) = -U(2,:,:);
    % super_guh(8,:,:) = U(1,:,:).^2 - U(2,:,:).^2;
    % super_guh(9,:,:) = U(1,:,:).*U(2,:,:);
    % 
    m_eq(1,:,:)=Rho;
    m_eq(2,:,:)=Rho.*(-2+3*sum(U.*U,1));
    m_eq(3,:,:)=Rho.*(1-3*sum(U.*U,1));
    m_eq(4,:,:)=Rho.*U(1,:,:);
    m_eq(5,:,:)=-Rho.*U(1,:,:);
    m_eq(6,:,:)=Rho.*U(2,:,:);
    m_eq(7,:,:)=-Rho.*U(2,:,:);
    m_eq(8,:,:)=Rho.*(U(1,:,:).*U(1,:,:)-U(2,:,:).*U(2,:,:));
    m_eq(9,:,:)=Rho.*(U(1,:,:).*U(2,:,:));
    % m_eq = Rho.*super_guh;
    % 
    % phi = -pagemtimes(pagemtimes(inv(M),S), (m-m_eq));
    % 
    % f = f_new + phi;
    %f=f_new-pagemtimes(pagemtimes(inv(M),S),(pagemtimes(M,f_new)-m_eq));
    
    f=f_new-pagemtimes(pagemtimes(inv(M),S),(pagemtimes(M,f_new)-m_eq));
    [guh, res_list(t)] = res(Rho_old, Rho, min_error);
    fprintf("Itt: %i     ||     Res: %.4e\n", t, res_list(t))
    %fprintf("Itt: %i\n", t)
end
total_time = toc;
%% Post-Processing / Visualization
clc;
load Ghia_Re100.mat
figure
quiver(flipud(squeeze(U(1,:,:))),flipud(squeeze(U(2,:,:))),10)
axis equal tight

figure
contourf(flipud(squeeze(Rho)),30)
axis equal tight

mid = N_x/2;
Vertical_Sample = U(1, :, mid)/U_top;
Horizontal_Sample = U(2, mid, :)/U_top;

%u2_Ghia = [1 0.48223 0.46120 0.45992 0.46036 0.33556 0.20087 0.08183 -0.03039 -0.07404 -0.22855 -0.33050 -0.40435 -0.43643 -0.42901 -0.41165 0.00000 ];
%v2_Ghia = [0.00000 -0.49774 -0.55069 -0.55408 -0.52876 -0.41442 -0.36214 -0.30018 0.00945 0.27280 0.28066 0.35368 0.42951 0.43648 0.43329 0.42447 0.00000];

figure
plot(squeeze(Vertical_Sample), flip((1:L)/L), "black", flip(u_Ghia), flip(y_Ghia))
title("Vertical Sample (U)")
xlabel("u")
ylabel("y")
figure
plot((1:L)/L, squeeze(Horizontal_Sample), "black", flip(x_Ghia), flip(v_Ghia));
title("Horizontal Sample (V)")
xlabel("x")
ylabel("v")

figure
u = flip(squeeze(U(1, :, :)));
v = squeeze(U(2, :, :));
[startX, startY] = meshgrid(1:50:N_x, 1:50:N_y);
verts = stream2(1:N_x,1:N_y,u,v,startX,startY);
streamline(verts)