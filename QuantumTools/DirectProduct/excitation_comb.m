function [excitations, indices] = excitation_comb(Nexcitations, Nsystems)
% The function computes the combinations of partition of several quantum
% excitations into several quantum systems. The function applies when the
% maximal allowed number of excitations in each system isn't smaller than the
% number of excitations. The algorithm is based on a trick which we learned
% in the statistical thermodynamics course of Avinoam Ben-Shaul in the
% tirgul.
% Input:
% Nexcitations: The number of excitations
% Nsystems: The number of systems
% Output:
% excitations: A matrix which represents the possible combinations of 
% partition of the excitations into the systems. The row index represents
% different combinations, and the column index represents the different
% systems.
% indices: A matrix which contains the indices of the vector representations
% of the combinations contained in the variable "excitations" for each of the
% individual systems. The structure is the same as that of "excitations".
% Note: The indices in the direct product space can be computed by the
% function kron_index, where the output variable indices in the present
% function is the second input variable in kron_index.
    % The number of partitions which divide the "balls" into separate
    % "baskets":
    Npartitions = Nsystems - 1;
    Ntotal = Nexcitations + Npartitions;
    % The vector 1:Ntotal represents the indices of the total elements, which
    % include the balls and the partitions. We compute the number of
    % combinations of choosing Npartitions partitions from the total
    % elements:
    partition_comb = nchoosek(1:Ntotal, Npartitions);
    % The number of combinations:
    Ncomb = size(partition_comb, 1);
    % We add fixed partitions at the beginning and the end, to define
    % the "baskets":
    partition_comb = [zeros(Ncomb, 1), partition_comb, ones(Ncomb, 1)*(Ntotal + 1)];
    indices = partition_comb(:, 2:(Nsystems + 1)) - partition_comb(:, 1:Nsystems);
    excitations = indices - 1;
end