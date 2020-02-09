%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author:       Pim van Hees
% Insitution:   Eindhoven University of Technology
% Department:   Mechanical Engineering
% Group:        Mechanics of Materials
% Subject:      4EM30 Multiscle Modelling for Polymer Mechanics
% Date:         08-02-2020
% Title:        Potential energy calculation for particle bonds
% Description:
%   Calculates the potential energy for a system of bonded particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Epot = calc_EpotBond(r,bond,k)
    % input:
    %   r:      N*dim matrix containing the positions of the N particles
    %   vond:   Nbond*3 matrix containing for Nbond bonds two particle
    %           numbers and the relaxed bond length
    %   k:      bond stiffness for all bonds
    %
    % output:
    %   Epot:   The potential energy of the system
    
    Epot = 0;
    for i = 1:size(bond,1)
        dist = r(bond(i,1),:) - r(bond(i,2),:);
        elong = norm(dist)-bond(i,3);
        Epot = Epot + k/2*elong^2;
    end
end