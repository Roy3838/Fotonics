# propagar.jl
using FFTW
using DSP

function meshgrid(x, y)
    X = repeat(x', length(y), 1)
    Y = repeat(y, 1, length(x))
    return X, Y
end

# function propagar(U_0, z, KXs, KYs, k_0)
#     F_U = fftshift(U_0,(1,2))
#     # xd
#     # kz = k_0 .- 0.5 .* (KXs.^2 .+ KYs.^2) ./ k_0
#     kt = sqrt.(KXs.^2 .+ KYs.^2)
#     kz = sqrt.(k_0^2 .- kt.^2)

#     F_Um = exp.(1im .* kz .* z) .* F_U
#     U = ifft(F_Um,(1,2))
#     return U
# end

function propagar(U_0, dz, KXs, KYs, k_0)
    # Calculate angular spectrum
    F_U = fftshift(fft(U_0, (1, 2)))
    
    # Calculate transverse wavevector magnitude
    kt = sqrt.(KXs.^2 .+ KYs.^2)
    
    # Use either the paraxial or non-paraxial propagator
    # Paraxial Propagator
    Prop = exp.(-1im * 0.5 * dz * kt.^2 / k_0)
    
    # Uncomment for non-paraxial propagation
    # kz = sqrt.(k_0^2 - kt.^2)
    # Prop = exp.(1im * dz * kz)
    
    # Propagation in frequency space
    F_U .*= Prop
    
    # Retrieve field in spatial domain
    U = ifft(F_U, (1, 2))
    
    return U
end





export meshgrid, propagar, conv_prop
