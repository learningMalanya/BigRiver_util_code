function transform_pheno_to_gemma(inputfile::Array{Float64, 2}, id::Int64, outputfile::AbstractString)
    pheno_y = inputfile[:, id];
    open(outputfile, "w") do io
        writedlm(io, pheno_y)
    end
end

function transform_geno_to_gemma(marker_values::Array{Float64, 2}, marker_names::Array{String, 1}, alleles::Array{String, 1}, outputfile::AbstractString)
    minor_allele = fill(alleles[1], length(marker_names), 1);
    major_allele = fill(alleles[2], length(marker_names), 1);
    output = hcat(hcat(marker_names, minor_allele, major_allele), transpose(marker_values))
    writeToFile(output, outputfile)
    return output
end

function gemmaWrapper(pheno_filename::String, geno_filename::String,
                      kinship_filename::String, output_filename::String)

    run(`$gemma -g $geno_filename -p $pheno_filename -k $kinship_filename -lmm 2 -lmax 1000000 -o $output_filename`)

end

function run_gemma_bxd(pheno::Array{Float64, 2}, geno::Array{Float64, 2}, kinship::Array{Float64, 2},
                       alleles::Array{String, 1}, marker_names::Array{String, 1},
                       gemma_path::String)

    (n, m) = size(pheno);
    p = size(geno, 2);

    gemma = gemma_path;

    L = zeros(p, m);

    writedlm("data/GEMMA_data/bxd_kinship.txt", kinship, '\t');
    transform_geno_to_gemma(geno, marker_names, alleles, "data/GEMMA_data/bxd_geno.txt");

    for i in 1:m

        transform_pheno_to_gemma(pheno, i, "data/GEMMA_data/bxd_pheno.txt");
        gemmaWrapper("data/GEMMA_data/bxd_pheno.txt", "data/GEMMA_data/bxd_geno.txt", "data/GEMMA_data/bxd_kinship.txt", "results_univariate_LMM");
        gemma_results = readdlm("output/results_univariate_LMM.assoc.txt", '\t');
        gemma_pvals = gemma_results[2:end, end] |> x -> Array{Float64}(x);
        gemma_lods = p2lod.(gemma_pvals, 1);
        L[:, i] = gemma_lods;

    end

    return L

end

function run_gemma(pheno::Array{Float64, 2}, geno::Array{Float64, 2}, kinship::Array{Float64, 2},
                       alleles::Array{String, 1}, marker_names::Array{String, 1},
                       pheno_filename::String, geno_filename::String, kinship_filename::String, output_filename::String, 
                       gemma_path::String)

    (n, m) = size(pheno);
    p = size(geno, 2);

    gemma = gemma_path;

    L = zeros(p, m);

    writedlm(kinship_filename, kinship, '\t');
    transform_geno_to_gemma(geno, marker_names, alleles, geno_filename);

    for i in 1:m

        transform_pheno_to_gemma(pheno, i, "$pheno_filename");
        gemmaWrapper("$pheno_filename", "$geno_filename", "$kinship_filename", "$output_filename");
        gemma_results = readdlm(joinpath("output/", "$output_filename", ".assoc.txt"), '\t');
        gemma_pvals = gemma_results[2:end, end] |> x -> Array{Float64}(x);
        gemma_lods = p2lod.(gemma_pvals, 1);
        L[:, i] = gemma_lods;

    end

    return L

end