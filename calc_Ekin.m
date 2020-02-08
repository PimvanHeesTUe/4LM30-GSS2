%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author:       Pim van Hees
% Insitution:   Eindhoven University of Technology
% Department:   Mechanical Engineering
% Group:        Mechanics of Materials
% Subject:      4EM30 Multiscle Modelling for Polymer Mechanics
% Date:         08-02-2020
% Title:        Kinetic energy calculation
% Description:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ekin = calc_Ekin(v,m)
    Ekin = 0;
    for i = 1:size(v,1)
        Ekin = Ekin + m/2*v(i,:)*v(i,:)';
    end
end