%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author:       Pim van Hees
% Insitution:   Eindhoven University of Technology
% Department:   Mechanical Engineering
% Group:        Mechanics of Materials
% Subject:      4EM30 Multiscle Modelling for Polymer Mechanics
%               Guided selfstudy 2
% Date:         07-02-2020
% Title:        Simulation of single polymer chain in 3D
% Description:     
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

%% initialisation
n = ceil(t_end/dt)+1; % amount of time steps

% plotted variables
pos = zeros(N,dim,n);           % position vectors
vel = zeros(N,dim,n);           % velocity vectors
Ekin = zeros(n,1);              % Kinetic energy
Epot = zeros(n,1);              % Potential energy
Etot = zeros(n,1);              % Total energy

% bonds between particles
bond = zeros(N-1,3);
for i = 1:N-1
    bond(i,:)=[i,i+1,l0];
end

% initial conditions
pos(:,1,1) = linspace(0,(N-1)*l0,N);
vel(:,:,1) = randn(N,dim)*velrms;
Ekin(1) = calc_Ekin(vel(:,:,1),m);
Epot(1) = calc_EpotBond(pos(:,:,1),bond,k);
Etot(1) = Ekin(1);

% force
Fnew = zeros(N,dim);

%% time looping
figure(1)
for i = 1:n-1
    % old force
    Fold = Fnew;
    % update position
    pos(:,:,i+1) = VelVerletPos(pos(:,:,i),vel(:,:,i),Fold,m,dt);
    % new force
    Fnew = forceall(pos(:,:,i+1),bond,k);
    % update velocity
    vel(:,:,i+1) = VelVerletVel(vel(:,:,i),Fold,Fnew,m,dt);
    
    % Energy
    Ekin(i+1) = calc_Ekin(vel(:,:,i+1),m);
    Epot(i+1) = calc_EpotBond(pos(:,:,i+1),bond,k);
    Etot(i+1) = Ekin(i+1) + Epot(i+1);
    
    % plotting
    figure(1)
    plot3(pos(:,1,i+1),pos(:,2,i+1),pos(:,3,i+1),'--ro')
    xlim([-5,15])
    ylim([-5,5])
    zlim([-5,5])
end

t=0:dt:t_end;
figure(2)
hold on
plot(t,Ekin,':r')
plot(t,Epot,':b')
plot(t,Etot,'-k')
legend({"Kinetic energy","Potential energy","Total energy"})

