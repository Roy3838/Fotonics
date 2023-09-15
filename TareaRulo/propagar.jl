using FFTW
using DSP

"""
    meshgrid(x, y)

Create 2D grid arrays similar to MATLAB's `meshgrid`.
"""
function meshgrid(x, y)
    X = repeat(x', length(y), 1)
    Y = repeat(y, 1, length(x))
    return X, Y
end

"""
    propagar(U_0, z, dz, KXs, KYs, k_0)

Propagate an input field `U_0` over a distance `z` using the Fourier transform method in steps of `dz`.

Parameters:
- `U_0`: Initial field to propagate.
- `z`: Total propagation distance.
- `dz`: Propagation step.
- `KXs`, `KYs`: Meshgrids of the spatial frequency in x and y directions.
- `k_0`: Wavenumber (2π/λ).

Returns:
- An array of propagated fields for each step `dz`.
"""
function propagar(U_0, z, dz, KXs, KYs, k_0)
    
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

    # Total number of propagation steps
    n_steps = Int(z/dz)
    
    # Array to hold all propagated fields
    propagated_fields = Array{Complex{Float64}, 3}(undef, size(U_0, 1), size(U_0, 2), n_steps+1)
    
    # Initial field
    propagated_fields[:, :, 1] = U_0
    # Loop through each propagation step
    for i in 1:n_steps
            
        # Propagation in frequency space
        F_U .*= Prop
        
        # Retrieve field in spatial domain
        U = ifft(F_U, (1, 2))
        
        # Store the propagated field in the array
        propagated_fields[:, :, i+1] = U
    end

return propagated_fields
end
