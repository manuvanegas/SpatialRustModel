
function in_farm(coord::CartesianIndex, side::Int)::Bool
    @inbounds for d in 1:2
        1 <= coord[d] <= side || return false
    end
    return true
end

function create_shade_map(farm_map::Matrix{Int}, shade_r::Int, side::Int, common_map::Symbol)

    shade_map = zeros(size(farm_map))
    
    if common_map != :fullsun
        maxdist = (2 * shade_r)^2
        crown = Iterators.filter(c -> c != CartesianIndex(0,0), CartesianIndices((-shade_r:shade_r, -shade_r:shade_r)))

        for coord in CartesianIndices(farm_map)
            if farm_map[coord] == 2
                shade_map[coord] = 1.0
            else
                neighs = (n for n in Iterators.filter(r -> in_farm(coord + r, side),
                    crown) if farm_map[coord + n] == 2
                )
                if !isempty(neighs)
                    shade_map[coord] = 1.0 - minimum(eucdist(n) for n in neighs) / maxdist
                end
            end
        end

    end

    return shade_map
end

eucdist(d::CartesianIndex{2})::Float64 = d[1]^2 + d[2]^2
