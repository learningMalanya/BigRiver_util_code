# This file contains the utility code to manipulate on the kinship matrix

## Code to construct kinship matrix for individual samples from a smaller number of strains, 
##  given the kinship matrix for strains and the samples-per-strain structure.

### helper function for calculating the kinship matrix for individual samples
function calcRepeats(x::Array{<:Any, 1})
    
    counting_dict = Dict{typeof(x[1]), Int64}();
    
    for i in 1:length(x)
        curr_key = x[i];
        if haskey(counting_dict, curr_key) # if current key already exists
            counting_dict[curr_key] += 1;
        else # if current key does not present, add it as a new key
            counting_dict[curr_key] = 1;
        end
    end
    
    return counting_dict
    
end

"""
    calcIndKinship_from_StrainKinship()

Construct kinship matrix for individual samples from strains given kinship for strains and the sample structure, 
    by duplicating elements of the original kinship matrix with respect to the individual sample structure.

# Arguments

- kinship_strains = m-by-m kinship matrix of m strains
- strain_info_about_samples = a vector of information about which strain each sample is from 

# Value

- kinship_ind_from_strains = the n-by-n kinship matrix of n individual samples

"""
function calcIndKinship_from_StrainKinship(kinship_strains::Array{Float64, 2}, 
                                           strain_info_about_samples::Array{<:Any, 1})
    
    counting_dict = calcRepeats(strain_info_about_samples);
    reps_each_strain = map(x -> counting_dict[x], unique(strain_info_about_samples));
    
    # Initialize the placeholder:
    kinship_ind_from_strains = ones(sum(reps_each_strain), sum(reps_each_strain));

    # process the first BXD strain...
    reps_strain_1 = reps_each_strain[1];
    k_11 = kinship_strains[1, 1];

    row_id = 0;
    col_id = 0;

    @views kinship_ind_from_strains[(row_id+1):(row_id+reps_strain_1), 
                             (row_id+1):(row_id+reps_strain_1)] .*= k_11;

    for j in 2:length(reps_each_strain)
        reps_strain_j = reps_each_strain[j];
        k_1j = kinship_strains[1, j];

        col_id += reps_each_strain[j-1];
        @views kinship_ind_from_strains[(row_id+1):(row_id+reps_strain_1), 
                                 (col_id+1):(col_id+reps_strain_j)] .*= k_1j;

        # process the off-diagonal block conveniently...
        kinship_ind_from_strains[(col_id+1):(col_id+reps_strain_j),
                                 (row_id+1):(row_id+reps_strain_1)] .*= k_1j;

    end

    # process for the second BXD strain and all after...
    for i in 2:length(reps_each_strain)

        reps_strain_i = reps_each_strain[i];
        k_ii = kinship_strains[i, i];

        row_id += reps_each_strain[i-1];
        @views kinship_ind_from_strains[(row_id+1):(row_id+reps_strain_i), 
                                 (row_id+1):(row_id+reps_strain_i)] .*= k_ii;

        col_id = row_id+reps_strain_i;

        for j in (i+1):length(reps_each_strain)

            reps_strain_j = reps_each_strain[j];
            k_ij = kinship_strains[i, j];

            @views kinship_ind_from_strains[(row_id+1):(row_id+reps_strain_i), 
                                     (col_id+1):(col_id+reps_strain_j)] .*= k_ij;

            # process the off-diagonal block conveniently...
            @views kinship_ind_from_strains[(col_id+1):(col_id+reps_strain_j),
                                     (row_id+1):(row_id+reps_strain_i)] .*= k_ij;

            col_id += reps_each_strain[j];
        end
    
    end
    
    return kinship_ind_from_strains
end

#####################################################################################################################
############################################## Alternative ##########################################################
#####################################################################################################################
function createRanges(values::Array{Int64, 1})
    
    m = length(values); # number of strains
    
    list_ranges = Array{UnitRange{Int64}, 1}(undef, m)
    startpt = 1
    
    for i in 1:length(values)
        endpt = startpt+values[i]-1;
        list_ranges[i] = startpt:endpt
        startpt = endpt+1
    end
    
    return list_ranges
    
end

function mapIds(id::Int64, ranges::Array{UnitRange{Int64}, 1})
    
    curr_group = 1;
    
    while !(id in ranges[curr_group])
        curr_group = curr_group+1;
    end
    
    return curr_group
    
end

function mapValues(idx::Int64, idy::Int64, ranges::Array{UnitRange{Int64}, 1}, K::Array{Float64, 2})
    
    f_idx = mapIds(idx, ranges);
    f_idy = mapIds(idy, ranges);
    
    return K[f_idx, f_idy];
    
end

function calcKinship2(kinship_strains::Array{Float64, 2}, 
    strain_info_about_samples::Array{<:Any, 1})

    counting_dict = calcRepeats(strain_info_about_samples);
    reps_each_strain = map(x -> counting_dict[x], unique(strain_info_about_samples));
    ranges = createRanges(reps_each_strain);

    n = sum(reps_each_strain); # total number of individual samples

    kinship_ind = ones(n, n);

    for i in 1:n
        for j in 1:n
            kinship_ind[i, j] = mapValues(i, j, ranges, kinship_strains);
        end
    end

    return kinship_ind;    

end