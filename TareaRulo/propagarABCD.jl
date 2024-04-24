using FFTW
using DSP

function meshgrid(x, y)
    X = repeat(x', length(y), 1)
    Y = repeat(y, 1, length(x))
    return X, Y
end

function propagar(U_0, z, dz, KXs, KYs, k_0, ABCDs=[], ABCD_zs=[])
    # Calculate angular spectrum
    F_U = fftshift(fft(U_0, (1, 2)))
    
    # Calculate transverse wavevector magnitude
    kt = sqrt.(KXs.^2 .+ KYs.^2)
    
    # Paraxial Propagator
    Prop = exp.(-1im * 0.5 * dz * kt.^2 / k_0)
    
    n_steps = Int(z/dz)
    propagated_fields = Array{Complex{Float64}, 3}(undef, size(U_0, 1), size(U_0, 2), n_steps+1)
    propagated_fields[:, :, 1] = U_0
    
    # Prepare a dictionary for ABCD matrix application at specific z positions
    ABCD_dict = Dict(ABCD_zs .=> ABCDs)

    for i in 1:n_steps
        current_z = i * dz
        
        # Propagation in frequency space
        F_U .*= Prop
        
        # Check for ABCD matrix application
        if current_z in keys(ABCD_dict)
            # Apply the ABCD matrix transformation
            U = apply_abcd(ifft(F_U, (1, 2)), ABCD_dict[current_z], KXs, KYs, k_0)
            # Transform back to frequency domain for further propagation
            F_U = fftshift(fft(U, (1, 2)))
        else
            # Continue in frequency space
            U = ifft(F_U, (1, 2))
        end
        
        # Store the propagated field in the array
        propagated_fields[:, :, i+1] = U
    end

    return propagated_fields
end

function apply_abcd(U, ABCD, KXs, KYs, k_0)
    # Extract ABCD matrix elements
    A, B, C, D = ABCD[1, 1], ABCD[1, 2], ABCD[2, 1], ABCD[2, 2]
    
    # Compute the wavevector magnitude squared (kt^2)
    kt_squared = KXs.^2 + KYs.^2
    
    # The phase modification is based on the ABCD matrix parameters B and C
    # Phase shift from C (curvature modification, related to lenses or curved mirrors)
    phase_shift_C = -1im * C * kt_squared / (2 * k_0)
    
    # Propagation distance phase adjustment from B (distance through which the beam propagates)
    phase_shift_B = -1im * B * k_0
    
    # Applying the ABCD transformation to the field's phase in the Fourier domain
    U_transformed = U .* exp.(phase_shift_C .+ phase_shift_B)
    
    return U_transformed
end
