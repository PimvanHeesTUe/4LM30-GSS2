%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author:       Pim van Hees
% Insitution:   Eindhoven University of Technology
% Department:   Mechanical Engineering
% Group:        Mechanics of Materials
% Subject:      4EM30 Multiscle Modelling for Polymer Mechanics
%               Guided selfstudy 3
% Date:         11-02-2020
% Title:        Simulation of single polymer chain in 3D
% Description:    
    % Simulation of polymer chain with boundary conditions
    % Gives animation of chain movemet and plots of end to end distance,
    % Energies, and end reaction force versus time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc;

%% Parameters
% particles
N  = 10;    % amount of particles
m  = 1 ;    % mass of single particle
l0 = 1 ;    % initial bond length
k  = 1 ;    % bond stiffness

% time-stepping
dt    = 0.01;   % time step
t_end = 10  ;   % length of simulation

% general
dim    = 3  ;   % 3D
velrms = 0.3;   % starting velocity root mean squared

% functionality
animation = false;      % if true, animation of polymer chain is shown
scheme = "velverlet";   % Which integration scheme is used: "velverlet" (Velocity-Verlet) or "euler"

% bc
fixed = [1,N]; % Numbers of the particles are fixed
Fx = 0.5;      % Force in x-direction on last particle

% random
% rng(1000);    % Random number generator seed. Used for comparing results

%% initialisation
n = ceil(t_end/dt)+1; % amount of time steps

% plotted variables
pos = zeros(N,dim,n);           % position vectors
vel = zeros(N,dim,n);           % velocity vectors
Ekin = zeros(n,1);              % Kinetic energy
Epot = zeros(n,1);              % Potential energy
Etot = zeros(n,1);              % Total energy
end2end = zeros(n,1);           % end to end distance of the chain
Fend = zeros(n,1);              % Reaction force on the last particle
Work = zeros(n,1);              % Work exerted by the force on the last particle
Ekin_part = zeros(N,1);         % Kinetic energy per particle integrated over time

% bonds between particles
bond = zeros(N-1,3);
for i = 1:N-1
    bond(i,:)=[i,i+1,l0];
end

% boundary conitions
bc_pos = false(N,dim);
bc_F = zeros(N,dim);
if ~isempty(fixed)
    bc_pos(fixed,:) = true(length(fixed),3);
end
bc_F(N,:) = [Fx,0,0];

% initial conditions
pos(:,1,1) = linspace(0,(N-1)*l0,N);
vel(:,:,1) = randn(N,dim)*velrms;
Ekin(1) = calc_Ekin(vel(:,:,1),m);
Epot(1) = calc_EpotBond(pos(:,:,1),bond,k);
Etot(1) = Ekin(1);
end2end(1) = norm(pos(1,:,1)-pos(end,:,1));
vel(:,:,1) = vel(:,:,1)-vel(:,:,1).*bc_pos;

% force
Fnew = bc_F;

%% time looping
% show the moving chain during the simulation
if animation
    figure(1)
end

for i = 1:n-1
    % old force
    Fold = Fnew;
    % update position
    pos(:,:,i+1) = VelVerletPos(pos(:,:,i),vel(:,:,i),Fold,m,dt);
    pos(:,:,i+1) = pos(:,:,i+1)-bc_pos.*(pos(:,:,i+1)-pos(:,:,i));
    % new force
    Fnew = forceall(pos(:,:,i+1),bond,k)+bc_F-bc_F.*bc_pos;
    % update velocity
    if strcmp(scheme,"velverlet")
        vel(:,:,i+1) = VelVerletVel(vel(:,:,i),Fold,Fnew,m,dt);
    elseif strcmp(scheme,"euler")
        vel(:,:,i+1) = VelVerletVel(vel(:,:,i),Fold,Fold,m,dt);
        % Velocity-Verlet scheme with fnew=fold for velocity computation is
        % the same as Euler scheme
    else
        fprintf("No such scheme!")
        return
    end
    
    vel(:,:,i+1) = vel(:,:,i+1)-bc_pos.*vel(:,:,i+1); % set velocities at fixed nodes to 0;
    
    % Energy
    Ekin(i+1) = calc_Ekin(vel(:,:,i+1),m);
    Epot(i+1) = calc_EpotBond(pos(:,:,i+1),bond,k);
    Work(i+1) = Fx*(pos(N,1,i+1)-pos(N,1,1));
    Etot(i+1) = Ekin(i+1) + Epot(i+1) - Work(i+1);
    for j = 1:N
        Ekin_part(j) = Ekin_part(j) + 1/2*m*vel(j,:,i+1)*vel(j,:,i+1)'*dt;
    end
    
    % End to end distance
    end2end(i+1) = norm(pos(1,:,i+1)-pos(end,:,i+1));
    
    % Force on the end of the chain
    Fend(i+1) = norm(Fnew(N,:));
    
    % plot the movement of the chain
    if animation
        figure(1)
        plot3(pos(:,1,i+1),pos(:,2,i+1),pos(:,3,i+1),'--ro') %#ok<*UNRCH>
        xlim([-5,15])
        ylim([-5,5])
        zlim([-5,5])
    end
end

t=0:dt:t_end;
figure(2)
hold on
plot(t,Ekin,'--r')
plot(t,Epot,'--b')
plot(t,Work,'--g')
plot(t,Etot,'-k')
legend({'Kinetic energy','Potential energy','Work','Total energy'})
title('Energy of the system')
xlabel('time')
ylabel('Engergy')

figure(3)
plot(t,end2end,'r')
xlabel('time')
ylabel('distance')
title('End to end distance of the chain')

figure(4)
plot(t,Fend,'r')
xlabel('time')
ylabel('force')
title('Normalized reaction force on the last particle')

for j = 1:N
    fprintf('Time averaged kinetic energy of particle %3d is %10.3e \n',j,Ekin_part(j)/t_end);
end