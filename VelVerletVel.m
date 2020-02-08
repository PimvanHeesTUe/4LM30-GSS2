%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author:       Pim van Hees
% Insitution:   Eindhoven University of Technology
% Department:   Mechanical Engineering
% Group:        Mechanics of Materials
% Subject:      4EM30 Multiscle Modelling for Polymer Mechanics
% Date:         08-02-2020
% Title:        Velocity-Verlet algorithm for velocity
% Description:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vnew] = VelVerletVel(vold,Fold,Fnew,m,dt)   
    vnew = vold + (Fnew + Fold)/(2*m)*dt;
end