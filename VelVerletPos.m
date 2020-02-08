%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author:       Pim van Hees
% Insitution:   Eindhoven University of Technology
% Department:   Mechanical Engineering
% Group:        Mechanics of Materials
% Subject:      4EM30 Multiscle Modelling for Polymer Mechanics
% Date:         08-02-2020
% Title:        Velocity-Verlet algorithm for position
% Description:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnew] = VelVerletPos(rold,vold,Fold,m,dt)   
    rnew = rold + vold*dt + Fold/(2*m)*dt^2;
end
