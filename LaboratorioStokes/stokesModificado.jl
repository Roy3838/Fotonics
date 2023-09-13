using Plots
using Images
using ImageFiltering
using Luxor

# Noise suppression function
function noise_supp(img::Array{Float64,2}, sze::Int)::Array{Float64,2}
    kern1 = 1/(2*sze+1) * ones(1,2*sze+1)
    kernf = kernelfactors((kern1, kern1))
    return imfilter(img, kernf)
end

# Crop function
function crop(img::Array{Float64,2}, x1::Int, x2::Int, y1::Int, y2::Int)::Array{Float64,2}
    return img[x1:x2, y1:y2]
end

# Path to images
roy_path = "C:/Users/maxga/Downloads/Photonics/LaboratorioStokes/ImagesCut/"
IM_PATH = roy_path

# Load images
Ix = Float64.(Gray.(load(IM_PATH*"1_1Cut.jpg")))
Iy = Float64.(Gray.(load(IM_PATH*"1_2Cut.jpg")))
Id = Float64.(Gray.(load(IM_PATH*"1_3Cut.jpg")))
Ia = Float64.(Gray.(load(IM_PATH*"1_4Cut.jpg")))
Ir = Float64.(Gray.(load(IM_PATH*"1_5Cut.jpg")))
Il = Float64.(Gray.(load(IM_PATH*"1_6Cut.jpg")))

# Compute Stokes parameters
S0 = Iy + Ix
S1 = Ix - Iy
S2 = Ir - Il
S3 = Id - Ia

# Crop the images
yc = 1:300
xc = 1:400
S0c = crop(S0, yc.start, yc.stop, xc.start, xc.stop)
S1c = crop(S1, yc.start, yc.stop, xc.start, xc.stop)
S2c = crop(S2, yc.start, yc.stop, xc.start, xc.stop)
S3c = crop(S3, yc.start, yc.stop, xc.start, xc.stop)

# Create the main plot with the image
I_0 = crop(Ix, yc.start, yc.stop, xc.start, xc.stop)
p = heatmap(I_0, color=:grays, axis=false)

# Overlay ellipses
N = 8
heightc = length(yc)
widthc = length(xc)
elp = palette([:red, :blue], 3)
for ii in 1:N:heightc
    for jj in 1:N:widthc
        S0_p = S0c[ii,jj]
        S1_p = S1c[ii,jj]
        S2_p = S2c[ii,jj]
        S3_p = S3c[ii,jj]

        #clr = Int((sign(S2_p))+2)
        #χ = 0.5 * atan(S3_p, S1_p)
        #b = sqrt((1 - atan((S1_p^2 + S3_p^2), S2_p^2) / (pi/2)) / 2)
        #a = sqrt(1 - b^2)
        clr = Int(sign(S3_p)+2)
        χ = -(atan.(S2_p,S1_p))./2 
        b = (S0_p .- hypot.(S1_p, S2_p))/2
        a = (S0_p .+ hypot.(S1_p, S2_p))/2 

        # Ellipse coordinates
        theta = range(0, stop=2*pi, length=100)
        x = a * cos.(theta) * sin(χ) - b * sin.(theta) * cos(χ)
        y = a * cos.(theta) * cos(χ) + b * sin.(theta) * sin(χ)

        # Scale and position
        xs = (N/2) .* x .+ jj
        ys = (N/2) .* y .+ ii

        plot!(xs, ys, color=elp[clr], lw=2, label=false)
    end
end

display(p)
