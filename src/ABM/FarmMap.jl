
function create_farm_map(map_side::Int = 100, row_d::Int = 2, plant_d::Int = 1, shade_d::Int = 6,
    shade_pattern::Symbol = :regular, barrier_rows::Int = 2, barriers::NTuple{2, Int} = (1, 0))::Array{Int,2}
    
    side = map_side
    # starting coffees
    farm_map = create_fullsun_farm_map(side, row_d, plant_d)

    # add shades
    if shade_d != 0 && (shade_d != side)
        if shade_pattern == :regular
            add_regshades!(farm_map, side, shade_d)
        else
            nshades = round(Int, (side / shade_d)^2)
            @inbounds farm_map[sample(1:side^2, nshades, replace=false)] .= 2
        end
    end

    midbarrs = @inbounds barriers[1] # internal barriers
    if midbarrs > 0
        barrow = barr_places(side, midbarrs, barrier_rows)
        barrcol = copy(barrow)

        if midbarrs == 1
            row_d == 3 && (barrcol .-= barrier_rows)
            row_d == 4 && barrier_rows == 2 && (barrcol .-= 1)
        elseif midbarrs > 1
            if barrier_rows == 1
                if plant_d == 2
                    # if suggested placement is odd, change it to even to place
                    # shades between coffee plants instead of replacing them
                    @inbounds barrow[isodd.(barrow)] .+= 1
                end
                if row_d > 1 # farm rows are array cols here
                    @inbounds barrcol[findall(x -> (x % row_d == 1), barrcol)] .+= 1
                end
            else
                if row_d > 2
                    conflicting = findall(x -> (x % row_d == 1), barrcol)
                    for cn in conflicting
                        # determine if cn is in 1st or 2nd half of placements:
                        # 1st half -> initial x's (equal to the result when barrier_rows is 1) 
                        # 2nd half -> extra x's because barrier_rows is 2
                        # eg, barr_places(100,2,1) = [33,66]; barr_places(100,2,2) = [33,66,34,67]
                        @inbounds if cn <= (length(barrcol) / 2)
                            @inbounds barrcol[[cn, (cn + midbarrs)]] .+= 1
                        else
                            @inbounds barrcol[[cn, (cn - midbarrs)]] .-= 1
                        end
                    end
                end
            end
        end

        @inbounds farm_map[barrow, :] .= 2
        @inbounds farm_map[:, barrcol] .= 2
    end

    if barriers[2] == 1
        @inbounds farm_map[[1, side], :] .= 2
        @inbounds farm_map[:, [1, side]] .= 2
    end

    return farm_map
end

function create_fullsun_farm_map(side::Int, row_d::Int, plant_d::Int)::Array{Int,2}
    farm_map = zeros(Int, side, side)
    for r in 1:row_d:side # these are farm "rows", but in the array they are columns (heatmap shows them as rows again)
        for p in 1:plant_d:side
            @inbounds farm_map[p, r] = 1
        end
    end
    return farm_map
end

function create_regshaded_farm_map(side::Int, row_d::Int, plant_d::Int, shade_d::Int)::Array{Int,2}
    farm_map = create_fullsun_farm_map(side, row_d, plant_d)
    for si in 1:shade_d:side
        for sj in 1:shade_d:side
            @inbounds farm_map[sj, si] = 2
        end
    end
    return farm_map
end

function add_regshades!(farm_map::Matrix{Int}, side::Int, shade_d::Int)::Array{Int,2}
    for si in 1:shade_d:side
        for sj in 1:shade_d:side
            @inbounds farm_map[sj, si] = 2
        end
    end
    return farm_map
end

## Helper
function barr_places(side::Int, num::Int, singledouble::Int)::Vector{Int}
    if num == 1
        placements = [fld(side, (num + 1))]
    elseif num ==2
        placements = [33, 70]
    end
    if singledouble == 2
        placements = vcat(placements, (placements .+ 1))
    end
    return placements
end

function adjust_pos(placements::Vector{Int}, idx::Vector{Int})
    v = copy(placements)
    v[idx] .+= 1
    return v
end

