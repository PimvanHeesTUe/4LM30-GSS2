%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author:       Pim van Hees
% Insitution:   Eindhoven University of Technology
% Department:   Mechanical Engineering
% Group:        Mechanics of Materials
% Subject:      4EM30 Multiscle Modelling for Polymer Mechanics
% Date:         08-02-2020
% Title:        Potential energy calculation for particle bonds
% Description:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Epot = calc_EpotBond(r,bond,k)
    Epot = 0;
    for i = 1:size(bond,1)
        dist = r(bond(i,1),:) - r(bond(i,2),:);
        elong = norm(dist)-bond(i,3);
        Epot = Epot + k/2*elong^2;
    end
end