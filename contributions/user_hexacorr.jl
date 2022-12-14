# * ==================================================
# * Purpose: DSF transfrmation for hexagonal unit cell
# * Written: Chaebin Kim (Seoul National University)
# * Date: 09/20/2022
# * ==================================================

# ? hexa! : Correcting the momentum space with hexagonal unit cell

"""
    user_hexa_corr!(MySQWperp)

Transformation for hexagonal lattice

"""

function user_hexacorr!(MySQWperp)
    A = OffsetArrays.no_offset_view(MySQWperp.sfactor) # Copying the sfactor wihtout offsets
    Lx, Ly, Lz, T = size(A)

    New = zeros(round(Int, 3Lx / 2), Ly, Lz, T)
    for i = 1:Lx, j = 1:Ly
        y_pot = div(i + 2 * j, 2, RoundUp)
        New[y_pot, i, :, :] .= A[j, i, :, :]
    end

    L_min = -1 .* div.(size(New) .- 1, 2)
    New_OSA = OffsetArray(New, OffsetArrays.Origin(L_min...)) # Adding offsets

    return New_OSA
end