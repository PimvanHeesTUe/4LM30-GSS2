%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author:       Pim van Hees
% Insitution:   Eindhoven University of Technology
% Department:   Mechanical Engineering
% Group:        Mechanics of Materials
% Subject:      4EM30 Multiscle Modelling for Polymer Mechanics
% Date:         08-02-2020
% Title:        Velocity-Verlet algorithm for position
% Description:  
%   Calculates the position in the new time step of particles based on the
%   velocity and force in the current time step.
%   Based on Velocity-Verlet scheme
%   However position calculation in the Euler scheme is identical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnew] = VelVerletPos(rold,vold,Fold,m,dt)
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
    %   rnew    N*dim matrix containing particle positions on the new
    %           time step
    
    rnew = rold + vold*dt + Fold/(2*m)*dt^2;
end
