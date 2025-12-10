import Pkg
Pkg.instantiate()

using Revise

Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

using PeriodicAutoregressive
@info("""
This session is using PeriodicAutoregressive.jl with Revise.jl.
For more information visit https://timholy.github.io/Revise.jl/stable/.
""")
