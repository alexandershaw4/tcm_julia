# Set up model parameters and data structures

ns = 1 # num regions
np = 8 # num cells per region
nk = 7 # num states per cell

using DataStructures

M = OrderedDict{Char,Array{Float64}}()
P = OrderedDict{Char,Array{Float64}}()

# Fill out M struct
M['x'] = zeros(1,8,7)

# Fill out Parameter struct
P['A'] = ones(ns,ns) # connectivity
P['B'] = ones(ns,ns) # trial-spec connectivity
P['C'] = ones(ns,1)  # inputs
P['D'] = [0 0]       # delays [intr] [extr]
P['E'] = [0]         # background activity
P['G'] = [0]         # trial spec intrinsics (off)
P['H'] = zeros(ns,np,nk)

# State contributions to output
P['J'] = zeros(ns,np,nk)
P['J'][1,[1 2 4 6],1] = (log(Diagonal([.2,.8, .2, .2]))).diag'
P['L'] = [4]

# Other parameters
P['S']  = [0]
P['T']  = zeros(ns,4)
P['V'] = [0 0 0 0 0 0 0 0]
