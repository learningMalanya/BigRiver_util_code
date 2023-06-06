## To be moved...
function plotLL(y::Array{Float64, 2}, G::Array{Float64, 2}, covar::Array{Float64, 2}, 
    K::Array{Float64, 2}, 
    h2_grid::Array{Float64, 1}, markerID::Int64;
    x_lims::Array{Float64, 1}, y_lims::Array{Float64, 1},
    prior::Array{Float64, 1} = [0.0, 0.0])


ell_null = zeros(length(h2_grid));
ell_alt = zeros(length(h2_grid));


num_of_covar = size(covar, 2);
(y0, X0, lambda0) = transform_rotation(y, [covar G], K; addIntercept = false);

for k in 1:length(h2_grid)
curr_h2 = h2_grid[k];
output = getLL(y0, X0, lambda0, num_of_covar, markerID, curr_h2; prior = prior);
ell_null[k] = output.ll_null;
ell_alt[k] = output.ll_markerID;
end

opt_ll_null = findmax(ell_null)[1];
opt_h2_null = h2_grid[findmax(ell_null)[2]];
opt_ll_alt = findmax(ell_alt)[1];
opt_h2_alt = h2_grid[findmax(ell_alt)[2]];

p = plot(h2_grid, ell_null, xlabel = "h2", ylabel = "loglik", label = "Null", color = "blue", legend=:bottomleft)
scatter!(p, [opt_h2_null], [opt_ll_null], label = "maxLL_null", color = "blue")
plot!(p, h2_grid, ell_alt, xlabel = "h2", ylabel = "loglik", label = ("Alt_$markerID"), color = "red")
scatter!(p, [opt_h2_alt], [opt_ll_alt], label = "maxLL_alt", color = "red");

plot!(p, ones(2).*opt_h2_null, [y_lims[1]-0.05, opt_ll_null], color = "blue", style = :dash, label = "")
plot!(p, [x_lims[1]-0.05, opt_h2_null], ones(2).*opt_ll_null, color = "blue", style = :dash, label = "")
annotate!(p, opt_h2_null, y_lims[1]-1.75, 
  text("$opt_h2_null", :blue, :below, 8))

plot!(p, ones(2).*opt_h2_alt, [y_lims[1]-0.05, opt_ll_alt], color = "red", style = :dash, label = "")
plot!(p, [x_lims[1]-0.05, opt_h2_alt], ones(2).*opt_ll_alt, color = "red", style = :dash, label = "")
annotate!(p, opt_h2_alt, y_lims[1]-1.75, 
  text("$opt_h2_alt", :red, :below, 8))



xlims!(p, (x_lims[1]-0.05, x_lims[2]+0.05))
ylims!(p, (y_lims[1]-0.05), y_lims[2]+0.05)
#= 
ylims!(p, (minimum([y_lims[1], minimum(ell_null), minimum(ell_alt)])-0.05, 
   maximum([y_lims[2], maximum(ell_null), maximum(ell_alt)])+0.05))

=#
return p

end