%%% ETH Zurich
%%% Finite Element Analysis in Biomedical Engineering Spring 2017
%%% Assignment 2: Solute transport through a cartilage tissue sample (2D)
%%% Explicit time marching, finite difference formulation
%%% Selver Pepic
%%% 29.04.2017. Version 1 (bulk of the code)
%%% 16.05.2017. Version 2 (additional plots added)
%%% 17.05.2017. Version 3 (minor changes, saving data to logfile added)

%% Physical parameters and simulation parameters
% Spatial domain
Lx = 10;        % sample length in x-direction (in mm)
Ly = 10;        % sample length in y-direction (in mm)
lx = 2;         % channel length in x-direction (in mm)
ly = 2;         % channel length in y-direction (in mm)

% Time domain
tmax = 48*60*60;    % total duration of simulation (in seconds)
tplot = [0 1 6 24 36 48]*60*60; % time points to be ploted

% Diffusion coefficients
Dx = 40E-6; % in mm^2/s
Dy = 20E-6; % in mm^2/s

% Spatial domain discretisation
n = 50;
nx = n;       % number of gridpoints in x direction
ny = n;       % number of gridpoints in y direction
dx = Lx/(nx-1); % grid spacing in x direction (in mm)
dy = Ly/(ny-1); % grid spacing in y direction (in mm)
x = 0:dx:Lx;    % array of x values (used when ploting)
y = 0:dy:Ly;    % array of y values (used when ploting)
% NOTE: number of gridpoints in x and y has been kept separate to have a
% more versatile code for later

% Time domain discretisation
nt = 2000;          % number of timesteps
dt = tmax/nt;       % length of time step
t = 0:dt:tmax-dt;      % array of time values (for plotting)
ntplot = round(tplot./tmax * nt);   % timesteps to be ploted
% correction of timestep values outside the valid range (0 < ntplot <= nt)
for it = 1:length(ntplot)
    if (ntplot(it)<=0)
        ntplot(it) = 1;
    end
    if (ntplot(it)>nt)
        ntplot(it) = nt;
    end
end
    
% Checking the stability criterion
% criterion is that in: C(x,y,t+1) = (1-2*bx-2*by) * C(x,y,t) + ..., 
    % the factor (1-2*bx-2*by) needs to be >= 0
% since dx=dy, with defining the anisotropy factor a=Dy/Dx=0.5 we have:
    % (1-2*bx-2*a*bx) >=0, which leads to:
    % bx <= 0.5/(1+a) = 0.5/1.5 = 1/3 = 0.333 OR by <= 0.166

bx = Dx*dt/(dx*dx);     % stability parameter
by = Dy*dt/(dy*dy);     % stability parameter
%bx = 0.334;             % manual input, used for exploring the stability
%by = bx*Dy/Dx;
if ((1-2*bx-2*by)<0)
    disp('WARNING: scheme unstable!');
    disp(['Multiplicative factor = ',num2str(1-2*bx-2*by)]);    
else
    disp('All fine, scheme should be stable.');
end

%% Initial and boundary conditions
% Center and inner boundary in "pixel" coordinates
xo = round(0.5*nx);     % center of the domain in "pixels"
yo = round(0.5*ny);     % center of the domain in "pixels"
wx = round(0.5*lx/Lx * nx);   % channel x half-width in "pixels"
wy = round(0.5*ly/Ly * ny);   % channel y half-width in "pixels"

% Initial conditions
    C = zeros(ny,nx,nt);   % concentration zero in whole domain, see below for inner boundary
    % NOTE: order of ny and nx is reversed on purpose since this keeps
    % orientations of physical and computational domains same

% Boundary conditions (at all time points, all in mol/mm^3)
    C(:,1,:)  = 0;      % left boundary
    C(:,nx,:) = 0;      % right boundary
    C(1,:,:)  = 0;      % bottom boundary
    C(ny,:,:) = 0;      % top boundary
    C(xo-wx:xo+wx, yo-wy:yo+wy,:) = 100;    % inner boundary (but also see below)

%    C(xo-wx:xo+wx, yo-wy:yo+wy,1:round(nt/2)) = 100;
%    C(xo-wx:xo+wx, yo-wy:yo+wy,round(nt(2):end) = 0;

% Plot of initial values (just for debbuging, not needed)
% figure
%    surf(x,y,C(:,:,1));
%    shading interp
%    title('Initial solute concentration (mol/mm^3)')
%    xlabel('x (mm)')
%    ylabel('y (mm)')
%    zlabel('Solute concentration (mol/mm^3)')
%    view(0,90);

%% FDM calculation
% Explicit time marching method, central differences
    % Cnew(0,0) = C(0,0) + Dx*dt/dx^2 * ( C(-1,0)+C(+1,0)-2*C(0,0) )
    %                    + Dy*dt/dy^2 * ( C(0,-1)+C(0,+1)-2*C(0,0) )
    % Cnew(0,0) = C(0,0) + bx * (C(-1,0)+C(+1,0)-2*C(0,0))
    %                    + by * (C(0,-1)+C(0,+1)-2*C(0,0))
    % Cnew(0,0) = (1-2bx-2by)*C(0,0) + bx*(C(-1,0)+C(+1,0))+ by*(C(0,-1)+C(0,+1));
    % Cnew(0,0) = (1-2b-2ab)*C(0,0) + b*(C(-1,0)+C(+1,0))+ ab*(C(0,-1)+C(0,+1));

% loop over internal nodes only, points in the channel are also included (for simplicity)
% since this would "diffuse away" the concentration in the channel, therefore
% the inner boundary BC is "reset" at each time step;
i = 2:nx-1;
j = 2:ny-1;
for k = 1:nt-1
    %if (k<=round(nt/2))
        C(xo-wx:xo+wx, yo-wy:yo+wy,k) = 100;
    %else
    %    C(xo-wx:xo+wx, yo-wy:yo+wy,k) = 0;
    %end
    C(j,i,k+1) = (1-2*bx-2*by) * C(j,i,k) +...
               bx * (C(j,i-1,k)+C(j,i+1,k)) + by * (C(j+1,i,k)+C(j-1,i,k));
end

%% Save data to a logfile
%   save( ['log_x',num2str(n),'_t',num2str(nt)] ,'C');
    
%% Ploting temperature profile at specified time points
for it = 1:length(ntplot)
    figure(it)
        surf(x,y,C(:,:,ntplot(it)));
        shading interp
        title({['Solute concentration (mol/mm^3) at t = ', num2str(tplot(it)/3600),' hours'],...
            ['dL = ',num2str(dx),' mm','; dT = ',num2str(dt),' s']})
        xlabel('x (mm)')
        ylabel('y (mm)')
        zlabel('Solute concentration (mol/mm^3)')
        view(0,90);
end

%% Ploting concentration profile along y = 1.4*yo (relative to the origin, not the center of channel)
yplot = 7;
ypix = round(yplot/Ly*ny);
figure
    plot(x,C(ypix,:,ntplot(5)));
    title({['Solute concentration (mol/mm^3) at t = ', num2str(tplot(5)/3600),' hours along y = ', num2str(yplot), ' mm'],...
            ['dL = ',num2str(dx),' mm','; dT = ',num2str(dt),' s']})
    xlabel('x (mm)')
    ylabel('Solute concentration (mol/mm^3)')
    
%% plot conc. vs time at non center point
plot(t,squeeze(C(25,25,:)))
hold on;
plot(t,squeeze(C(19,25,:)))
plot(t,squeeze(C(15,25,:)))
plot(t,squeeze(C(10,25,:)))
plot(t,squeeze(C(5,25,:)))
plot(t,squeeze(C(1,25,:)))