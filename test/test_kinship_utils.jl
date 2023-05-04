# Tests for utility functions for kinship matrices:

using Test
using LinearAlgebra

example_strains_num = [1, 3, 5, 7, 9, 11];
example_strains = "BXD" .* string.(example_strains_num);
example_samples_per_strain_num = vcat(example_strains_num[1:2], repeat(example_strains_num[3:6], inner = 3));
example_samples_per_strain = vcat(example_strains[1:2], repeat(example_strains[3:6], inner = 3));

rep_num = calcRepeats(example_samples_per_strain_num);
rep_string = calcRepeats(example_samples_per_strain);

example_K_strains = diagm(ones(6).*0.5).+ones(6, 6).*0.5


K_num = calcIndKinship_from_StrainKinship(example_K_strains, example_samples_per_strain_num)
K_string = calcIndKinship_from_StrainKinship(example_K_strains, example_samples_per_strain)


@test sum(abs.(K_num .- K_string)) <= 1e-10