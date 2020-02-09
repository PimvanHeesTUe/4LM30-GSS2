%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author:       Pim van Hees
% Insitution:   Eindhoven University of Technology
% Department:   Mechanical Engineering
% Group:        Mechanics of Materials
% Subject:      4EM30 Multiscle Modelling for Polymer Mechanics
% Date:         08-02-2020
% Title:        Velocity-Verlet algorithm for velocity
% Description: 
%   Calculates the velocity in the new time step of particles based on the
%   velocity and force in the current time step and the force in the new
%   time step.
%   Based on Velocity-Verlet scheme
%   However if for Fnew the value of Fold is used it becomes the Euler
%   scheme instead.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vnew] = VelVerletVel(vold,Fold,Fnew,m,dt)
    % input:
    %   rold:   N*dim matrix containing particle positions in current time step
    %   vold:   N*dim matrix containing particle velocities in current
    %           time step
    %   Fold:   N*dim matrix containing the point forces on the particles
    %           in the current time step
    %   m:      mass of all particles
    %   dt      time step
    %
    % output:
    %   vnew    N*dim matrix containing particle velocities on the new
    %           time step
    
    vnew = vold + (Fnew + Fold)/(2*m)*dt;
end