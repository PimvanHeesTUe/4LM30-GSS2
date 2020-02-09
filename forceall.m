%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author:       Pim van Hees
% Insitution:   Eindhoven University of Technology
% Department:   Mechanical Engineering
% Group:        Mechanics of Materials
% Subject:      4EM30 Multiscle Modelling for Polymer Mechanics
% Date:         07-02-2020
% Title:        Force computation for particles with bonds
% Description:
%   Calculates the forces on bonded particles by calculating the force per
%   bond. 
%   First calculates the distance and direction between bonded
%   particles and then uses this with the stiffness and relaxed bond length
%   to calculate the force vector on both particles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Fnew = forceall(pos,bond,k)
    % input:
    % pos: N*dim matrix of positions on current time step
    % bond: (N-1)*3 matrix containing (N-1) bonds with [particle 1,
    % particle 2, initial bond length]
    % k: stiffness
    %
    % output:
    % Fnew : N*dim matrix with current forces for all particles
    
    Fnew = zeros(size(pos,1),size(pos,2));
    
    for i = 1:size(bond,1)
        Fpair = forcepair(pos(bond(i,1),:),pos(bond(i,2),:),bond(i,3),k);
        Fnew(bond(i,2),:) = Fnew(bond(i,2),:) + Fpair;
        Fnew(bond(i,1),:) = Fnew(bond(i,1),:) - Fpair;
    end
end

function Fpair = forcepair(pos1,pos2,l0,k)
    dist = norm(pos2-pos1);    % current distance between the particles
    dir = (pos2-pos1)/dist;    % normalized direction
    Fabs = k*(l0 - dist);      % absolute force
    Fpair = Fabs * dir;     % force with direction for particle 2
end