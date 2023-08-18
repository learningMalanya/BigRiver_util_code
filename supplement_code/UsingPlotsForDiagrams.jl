# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# # Using Plots for Diagrams
# ___
# This notebooks explores how draw diagrams using Plots.

# + tags=[]
using Plots

# + [markdown] tags=[]
# ## Functions
# -

# ### Customized Shape

# + tags=[]
struct Matshape
    w
    h
    x
    y
end
# -

?Shape

# + tags=[]
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
rectangle(O::Matshape) = Shape(O.x .+ [0,O.w,O.w,0], O.y .+ [0,0,O.h,O.h])

# + tags=[]
# Returns 3 shapes: a matrix with one column and one row
matcolrow(A, B, C) = return rectangle(A), rectangle(B), rectangle(C)  
# -

# ### Plot functions

# + tags=[]
# M, C, R are respectively the 3 shapes: matrix, column, row
function plot_matcolrow!(M, C, R, my_color, mat_linewidth, mat_linecolor = my_color)
    # Plot matrix
    plot!(
        M,
        fillcolor = plot_color(my_color, 0.4), 
        linecolor = mat_linecolor,
        linewidth = mat_linewidth
    )

    # Plot matrix column
    plot!(
        C,
        fillcolor = plot_color(my_color, 0.3), 
        linecolor = mat_linecolor,
        linewidth = 1
    )

    # Plot matrix row
    my_plot = plot!(
        R,
        fillcolor = plot_color(my_color, 0.3), 
        linecolor = mat_linecolor,
        linewidth = 1
    )

    return my_plot
end
# -

# ## MLM Figure

# ### Draw a rectangle

# + tags=[]
# Matshape(width, height, x_coord, y_coord);
Y = Matshape(65, -60, 0, 0);
X = Matshape(45, -60, 70, 0);
L = Matshape(65, -90, 0, -65);
gap = X.x-Y.w;

# + tags=[]
rec_Y = rectangle(Y)
rec_X = rectangle(X)
rec_L = rectangle(L)
# -

# ### Customized attributes

# + tags=[]
my_linewidth = 2
# -

# ### Customized Colors

# + tags=[]
my_orange = RGB(241/255,163/255,64/255)
my_blue = RGB(103/255,169/255,207/255)
my_red = RGB(250/255,68/255,48/255)
my_black = RGB(0/255,0/255,0/255)
my_grey = RGB(200/255,200/255,200/255)
my_white = RGB(1, 1, 1)
my_green = RGB(0.76,1.0,0.76)
# -

# ### Plot

# + tags=[]
#################
# Plot Matrices #
#################
pInit = plot(
    legend = false,
    grid = false,
    xaxis = false,
    yaxis = false,
    right_margin = (3,:mm),
    left_margin = (-10,:mm),
    rec_Y,
    fillcolor = plot_color(my_orange, 0.4), 
    linecolor = my_orange,
    linewidth = my_linewidth
)

plot!(
    rec_X,
    fillcolor = plot_color(my_blue, 0.4), 
    linecolor = my_blue,
    linewidth = my_linewidth
)

plot!(
    rec_L,
    fillcolor = plot_color(my_green, 0.4), 
    linecolor = my_green,
    linewidth = my_linewidth
)

####################
# Plot Annotations #
####################
# font attributes
mat_font = "New Century Schoolbook Bold"
mat_font_sup = "New Century Schoolbook Bold Italic"
dim_font = "New Century Schoolbook Italic"

mat_font_size = 16
dim_font_size = 14
# my_font = "Bookman Demi"

# Matrix Names
annotate!([
        (Y.w/2, Y.h/2, ("Y", mat_font_size, :black, :center, mat_font)), 
        (X.w/2+X.x, X.h/2, ("X", mat_font_size, :black, :center, mat_font)),
        (L.w/2, L.h/2+L.y, ("L", mat_font_size, :black, :center, mat_font)),
])

# Descriptions
annotate!([
        (Y.w/2, Y.h/2-9, ("expression traits to be scanned", dim_font_size, :grey, :center, dim_font)), 
        (X.w/2+X.x, X.h/2-9, ("genome markers", dim_font_size, :grey, :center, dim_font)),
        (L.w/2, L.h/2+L.y-9, ("result matrix of LOD scores", dim_font_size, :grey, :center, dim_font)),
])

# Dimensions
annotate!([
        (Y.w+X.w+gap+4, Y.h/2, ("n", dim_font_size, :black, :center, dim_font)), 
        (Y.w+X.w-gap-8, X.h-5, ("p", dim_font_size, :black, :center, dim_font)),
        (L.w+3, Y.h+L.h/2-gap, ("p", dim_font_size, :black, :center, dim_font)),
        (L.w/2, Y.h+L.h-10, ("m", dim_font_size, :black, :center, dim_font)),
])

savefig("images/mlmdiagram.svg")
pInit
# -

# ## MLM Statin Figure

# ### Draw a rectangle

# + tags=[]
Y_box = Matshape(65,-60,0,0);
Y_col = Matshape(5,-60,55,0);
Y_row = Matshape(65,-5,0,-15);

X_box = Matshape(30,-60,70,0)
X_col = Matshape(5,-60,90,0);
X_row = Matshape(30,-5,70,-15);

L_box = Matshape(65,-35,0,-65)
L_col = Matshape(5,-35,55,-65);
L_row = Matshape(65,-5,0,-75);

rec_Y_box, rec_Y_col, rec_Y_row = matcolrow(Y_box, Y_col, Y_row);
rec_X_box, rec_X_col, rec_X_row = matcolrow(X_box, X_col, X_row);
rec_L_box, rec_L_col, rec_L_row = matcolrow(L_box, L_col, L_row);


# + tags=[]
#################
# Plot Matrices #
#################
pExample = plot(
    legend = false,
    grid = false,
    xaxis = false,
    yaxis = false,
    xlim = (-34,105),
    ylim = (-100,15),
    right_margin = (3,:mm),
    left_margin = (-10,:mm),
)

#################
# Plot Y matrix #
#################

plot_matcolrow!(rec_Y_box, rec_Y_col, rec_Y_row, my_orange, my_linewidth)


#################
# Plot X matrix #
#################

plot_matcolrow!(rec_X_box, rec_X_col, rec_X_row, my_blue, my_linewidth)


#################
# Plot L matrix #
#################
plot_matcolrow!(rec_L_box, rec_L_col, rec_L_row, my_blue, my_linewidth)

####################
# Plot Annotations #
####################
# font attributes
mat_font = "New Century Schoolbook Bold"
mat_font_sup = "New Century Schoolbook Bold Italic"
dim_font = "New Century Schoolbook Italic"

mat_font_size = 16
dim_font_size = 12
example_font_size = 10
# my_font = "Bookman Demi"

# Matrix Names
nntt_mat(M, C, R, my_text, my_font_size, my_color, my_alginmen, my_font) = ((M.x+C.x)/2, (M.y+M.h+R.y+R.h)/2, (my_text, my_font_size, my_color, my_alginmen, my_font))
nntt_mat_sup(M, C, R, my_text, my_font_size, my_color, my_alginmen, my_font) = ((M.x+C.x)/2+2.75, (M.y+M.h+R.y+R.h)/2+2.75, (my_text, my_font_size, my_color, my_alginmen, my_font))

annotate!([
        nntt_mat(Y_box, Y_col, Y_row, "Y", mat_font_size, :black, :center, mat_font),
        nntt_mat(X_box, X_col, X_row, "X", mat_font_size, :black, :center, mat_font),
        nntt_mat(Z_box, Z_col, Z_row, "Z", mat_font_size, :black, :center, mat_font),
])

annotate!([
        nntt_mat_sup(Z_box, Z_col, Z_row, "T", 9, :black, :center, mat_font_sup),
        nntt_mat_sup(B_box, B_col, B_row, "T", 9, :black, :center, mat_font_sup),
])



# Descriptions
annotate!([
        (Y_col.w/2+Y_col.x, 12, ("TG\n(12:0_12:0_18:2)", example_font_size, :grey, :center, dim_font)),
        (-7, Y_row.h/2+Y_row.y, ("sample 7", example_font_size, :grey, :right, dim_font)),
        (X_col.w/2+X_col.x, 9, ("Fish Oil", example_font_size, :grey, :center, dim_font)),
        (-7, Z_row.h/2+Z_row.y, ("Double Bonds = 2", example_font_size, :grey, :right, dim_font)),
])

# Dimensions
annotate!([
        (Y_box.w+X_box.w+gap+4, Y_box.h/2, ("n", dim_font_size, :black, :center, dim_font)), 
        (Y_box.w+X_box.w+gap+4, B_box.h/2+B_box.y, ("q", dim_font_size, :black, :center, dim_font)),
        (Z_box.w/2, Y_box.h+Z_box.h-gap-4, ("m", dim_font_size, :black, :center, dim_font))
])

# Horizontal Arrows
plot!(
    [-6, -1, Inf,-6,-1],
    [Y_row.h/2+Y_row.y, Y_row.h/2+Y_row.y, Inf, Z_row.h/2+Z_row.y,Z_row.h/2+Z_row.y],
    arrow=true,
    color=:grey,
    linewidth=1,
    label=""
)

plot!(
    [Y_col.w/2+Y_col.x, Y_col.w/2+Y_col.x, Inf, X_col.w/2+X_col.x, X_col.w/2+X_col.x],
    [6, 1, Inf,6,1],
    arrow=true,
    color=:grey,
    linewidth=1,
    label=""
)

savefig("images/mlmdiagramexample.svg")
pExample

# -

