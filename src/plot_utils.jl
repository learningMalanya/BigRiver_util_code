## To be moved...
function plotLL(y::Array{Float64, 2}, G::Array{Float64, 2}, covar::Array{Float64, 2}, 
    K::Array{Float64, 2}, 
    h2_grid::Array{Float64, 1};
    markerID::Int = 0,
    reml::Bool = false,
    prior::Array{Float64, 1} = [0.0, 0.0])


    ell_null = zeros(length(h2_grid));
    ell_alt = zeros(length(h2_grid));


    num_of_covar = size(covar, 2);
    (y0, X0, lambda0) = transform_rotation(y, [covar G], K; addIntercept = false);

    for k in 1:length(h2_grid)
        curr_h2 = h2_grid[k];
        output = getLL(y0, X0, lambda0, num_of_covar, markerID, curr_h2; reml = reml, prior = prior);
        ell_null[k] = output.ll_null;
        ell_alt[k] = output.ll_markerID;
    end

    opt_ll_null = findmax(ell_null)[1];
    opt_h2_null = h2_grid[findmax(ell_null)[2]];
    opt_ll_alt = findmax(ell_alt)[1];
    opt_h2_alt = h2_grid[findmax(ell_alt)[2]];

    if markerID > 0
        lb = min(findmin(ell_null)[1], findmin(ell_alt)[1]);
        ub = max(opt_ll_null, opt_ll_alt);
        else
        lb = findmin(ell_null)[1];
        ub = opt_ll_null;
    end

    p = plot(h2_grid, ell_null, xlabel = "h2", ylabel = "loglik", label = "Null", color = "blue", legend=:bottomleft)
    scatter!(p, [opt_h2_null], [opt_ll_null], label = "maxLL_null", color = "blue")

    if markerID > 0
        plot!(p, h2_grid, ell_alt, xlabel = "h2", ylabel = "loglik", label = ("Alt_$markerID"), color = "red")
        scatter!(p, [opt_h2_alt], [opt_ll_alt], label = "maxLL_alt", color = "red");
    end

    plot!(p, ones(2).*opt_h2_null, [lb-0.05, opt_ll_null], color = "blue", style = :dash, label = "")
    plot!(p, [-0.05, opt_h2_null], ones(2).*opt_ll_null, color = "blue", style = :dash, label = "")
    annotate!(p, opt_h2_null, lb-0.25,
    text("$opt_h2_null", :blue, :below, 8))

    if markerID > 0
        plot!(p, ones(2).*opt_h2_alt, [lb-0.05, opt_ll_alt], color = "red", style = :dash, label = "")
        plot!(p, [-0.05, opt_h2_alt], ones(2).*opt_ll_alt, color = "red", style = :dash, label = "")
        annotate!(p, opt_h2_alt, lb-0.25, text("$opt_h2_alt", :red, :below, 8))
    end


    xlims!(p, -0.05, 1.05)
    ylims!(p, lb-0.05, ub+0.05)

    return (fig = p, ell_null = ell_null)

end