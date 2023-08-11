function grow_shades!(current::Books, rate::Float64, max::Float64)
    current.ind_shade += rate * (1.0 - current.ind_shade * inv(max)) * current.ind_shade
end

