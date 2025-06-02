# Please also see performance notes at https://github.com/SebKrantz/OptimalTransportNetworks.jl
using DataFrames, CSV, LinearAlgebra, Statistics, Plots
using OptimalTransportNetworks
import HSL_jll
include("helpers.jl")

countries = unique(SubString.(readdir("data/grid_network/country/"), 1, 3))
failed = ["RUS", "KAZ"]

for c in ["KAZ"] # countries # ["ITA", "GRC"]  
    # if c in ["RUS", "KAZ"]
    #     continue
    # end      
    # Read Undirected Graph
    edges = CSV.read("data/grid_network/country/$(c)_edges.csv", DataFrame)
    # histogram(edges.time_efficiency, bins=100)
    # histogram(edges.cost_kmh, bins=100)

    # Read Nodes Data
    nodes = CSV.read("data/grid_network/country/$(c)_nodes.csv", DataFrame)
    nodes.pop_cell /= 1000 # Convert to thousands
    describe(nodes)

    n = maximum([maximum(edges.from), maximum(edges.to)])
    if nrow(nodes) > n 
        println("n is $(n) but nrow(nodes) is $(nrow(nodes))")
        nodes = nodes[1:n, :]
    end

    # Create Adjacency Matrix
    adj_matrix = falses(n, n)
    for i in 1:size(edges, 1)
        adj_matrix[edges.from[i], edges.to[i]] = adj_matrix[edges.to[i], edges.from[i]] = true
    end

    n_nodes = 1000 # nrow(nodes)

    # Create Infrastructure Matrix: Following Graff (2024) = average speed in km/h: length of route is accounted for in cost function
    infra_matrix = zeros(n, n)
    for i in 1:size(edges, 1)
        infra_matrix[edges.from[i], edges.to[i]] = infra_matrix[edges.to[i], edges.from[i]] = edges.time_efficiency[i]
    end
    println(describe(edges.time_efficiency))

    # Create Iceberg Trade Cost Matrix. Graff (2024): 0.1158826 * log(edges.distance[i] / 1.609)
    iceberg_matrix = zeros(n, n)
    for i in 1:size(edges, 1)
        iceberg_matrix[edges.from[i], edges.to[i]] = iceberg_matrix[edges.to[i], edges.from[i]] = 0.1158826 * log(edges.sp_distance[i])
    end
    iceberg_matrix[iceberg_matrix .< 0] .= 0
    extrema(iceberg_matrix)

    # Create Infrastructure Building Cost Matrix: Following Graff (2024)
    infra_building_matrix = zeros(n, n)
    for i in 1:size(edges, 1)
        infra_building_matrix[edges.from[i], edges.to[i]] = infra_building_matrix[edges.to[i], edges.from[i]] = edges.cost_kmh[i]
    end
    extrema(infra_building_matrix)

    # Basic characteristics of the economy
    population = nodes.pop_cell
    population += (population .== 0) * 1e-5

    # Productivity Matrix 
    # Check largest pcities
    productivity = zeros(n, maximum(nodes.product))
    for i in 1:n
        productivity[i, nodes.product[i]] = nodes.cell_GDPC_const_2017_PPP[i] * 1e4
    end
    extrema(productivity)
    all(sum(productivity .> 0, dims = 2) .== 1) # Check for Armington parameterization

    J, N = size(productivity)

    # Infrastructure Bounds
    min_mask = infra_matrix  # min.(infra_matrix, adj_matrix .* 10)  
    max_mask = max.(infra_matrix, adj_matrix .* 70) # Max time efficiency 
    # Cost of all existing infrastructure
    K_base = sum(infra_building_matrix .* infra_matrix)
    # Adjustments because of buffer
    K_ctry = 0.0
    for i in 1:size(edges, 1)
        if edges.is_buff[i]
            min_mask[edges.from[i], edges.to[i]] = min_mask[edges.to[i], edges.from[i]] = edges.time_efficiency[i]
            max_mask[edges.from[i], edges.to[i]] = max_mask[edges.to[i], edges.from[i]] = edges.time_efficiency[i]
        else
            K_ctry += 2 * edges.cost_kmh[i] * edges.time_efficiency[i]
        end
    end 
    # sum((infra_building_matrix .* infra_matrix)[.!nodes.is_buff, .!nodes.is_buff])
    # Now setting budget
    K = K_base + K_ctry * 0.1 # 10% increase

    if K < sum(infra_building_matrix .* min_mask)
        error("Infrastructure budget is too small for the network")
    end
    if K > sum(infra_building_matrix .* max_mask)
        error("Infrastructure budget is too large for the network")
    end

    # Parameters
    alpha = 0.4 # Spending share on traded goods in utility (curvature parameter)
    gamma = 0.1 # F&S: 0.1 ; TG: 0.946; Parameter governing intensity of congestion in transport
    beta = 0.13 # F&S: 0.13 ; TG: 1.2446 * gamma; Parameter governing returns to scale in infrastructure investment
    gamma = beta^2/gamma # IRS case
    sigma = 3.8; # Elasticity of substitution parameter (3.8 = Armington)
    a = 1 # F&S: 1; TG: 0.7; Returns to scale to labor in production function Zn * Ln^a
    rho = 0 # 2 # inequality aversion

    # Initialise geography
    param = init_parameters(annealing = false, labor_mobility = false, cross_good_congestion = n_nodes <= 200, duality = n_nodes > 200,
                            a = a, sigma = sigma, N = N, alpha = alpha, beta = beta, gamma = gamma, rho = rho, nu = 1,
                            K = K, tol = 1e-5, min_iter = 15, max_iter = 45, verbose = true)

    graph = create_graph(param, type = "custom", x = nodes.pwx, y = nodes.pwy, 
                         adjacency = adj_matrix, # omega = max.(.!nodes.is_buff, 0.1),
                         Lj = population, Zjn = productivity, Hj = population .* (1-alpha)) # TG: I normalise this because the general utility function has a (h_j/(1-alpha))^(1-alpha) thing with it

    graph[:delta_i] = infra_building_matrix;
    graph[:delta_tau] = iceberg_matrix;

    # Naming conventions: if IRS -> "cgc_irs" or "irs_na" without annealing; if alpha = 0.1 -> add "_alpha01"; if with_ports = false -> add "_noport"; if frictions -> add "_bc" (border cost)
    filename = "$(c)_10perc_fixed_irs_sigma38_rho0" # adjust if sigma != 1.5
    println("File extension: $filename")

    # Recommended to use coin HSL linear solvers. See README of OptimalTransportNetworks.jl and Ipopt.jl
    param[:optimizer_attr] = Dict(:hsllib => HSL_jll.libhsl_path, :linear_solver => "ma57") # Use ma86 for optimal performance on big machine 

    res_stat = res_opt = nothing
    try
        # Solve allocation from existing infrastructure
        @time res_stat = optimal_network(param, graph, I0 = infra_matrix, solve_allocation = true, verbose = true)

        # Solve Optimal Network (this can take long - up to 48h)
        @time res_opt = optimal_network(param, graph, I0 = infra_matrix, Il = min_mask, Iu = max_mask, verbose = false)
        # # Run Annealing Separately (this can take long - up to 100h)
        # @time res_opt, model, recover_allocation = optimal_network(param, graph, I0 = infra_matrix, Il = min_mask, Iu = max_mask, verbose = false, return_model = 2)
        # @time res_opt = annealing(param, graph, res_opt[:Ijk], final_model = model, recover_allocation = recover_allocation, allocation = res_opt, verbose = true)
    catch e
        append!(failed, [c])
        println(e)
        continue
    end
    # Welfare gains
    println("\nWelfare gain: ", round((res_opt[:welfare] / res_stat[:welfare] - 1) * 100, digits = 2), "%\n")

    # Check: should be 1
    K / sum(res_opt[:Ijk] .* graph[:delta_i])

    # Plot Network
    display(plot_graph(graph, res_opt[:Ijk], node_sizes = log.(res_opt[:Cj]), node_color = [:purple, :grey][nodes.is_buff .+ 1]))
    display(plot_graph(graph, res_opt[:Ijk] - infra_matrix, node_sizes = log.(res_opt[:Cj]), node_color = [:purple, :grey][nodes.is_buff .+ 1]))

    # Saving: Nodes
    res_nodes = deepcopy(nodes)
    res_nodes.uj_orig = vec(res_stat[:uj]) 
    res_nodes.Lj_orig = vec(res_stat[:Lj])
    res_nodes.Cj_orig = vec(res_stat[:Cj])
    res_nodes.Dj_orig = vec(res_stat[n_nodes > 200 ? :Cj : :Dj])
    res_nodes.PCj_orig = vec(res_stat[:PCj])
    # res_nodes.welfare = res_opt.welfare; # Same as: sum(res_opt.Lj .* res_opt.uj)
    res_nodes.uj = vec(res_opt[:uj])
    res_nodes.Lj = vec(res_opt[:Lj])
    res_nodes.Cj = vec(res_opt[:Cj])
    res_nodes.Dj = vec(res_opt[n_nodes > 200 ? :Cj : :Dj])
    res_nodes.PCj = vec(res_opt[:PCj])
    for n in 1:N
        res_nodes[!, Symbol("Lj_$(n)")] = res_opt[:Ljn][:,n]
        res_nodes[!, Symbol("Dj_$(n)")] = res_opt[n_nodes > 200 ? :Cjn : :Djn][:,n]
        res_nodes[!, Symbol("Yj_$(n)")] = res_opt[:Yjn][:,n]
        res_nodes[!, Symbol("Pj_$(n)")] = res_opt[:Pjn][:,n]
    end
    res_nodes |> CSV.write("results/grid_network/country/nodes_results_$(filename).csv")

    # Saving: Graph / Edges
    res_edges = deepcopy(edges)
    res_edges.Ijk_orig = res_to_vec(infra_matrix, edges)
    res_opt[:Ijk][nodes.is_buff, nodes.is_buff] .= infra_matrix[nodes.is_buff, nodes.is_buff]
    res_edges.Ijk = res_to_vec(res_opt[:Ijk], edges)
    for n in 1:N
        res_edges[!, Symbol("Qjk_$(n)")] = res_to_vec(res_opt[:Qjkn][:,:,n], edges)
    end
    res_edges |> CSV.write("results/grid_network/country/edges_results_$(filename).csv")
end




