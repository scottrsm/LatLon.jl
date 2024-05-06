using LatLon
import Pkg

Pkg.add("Documenter")
using Documenter

makedocs(
	sitename = "LatLon",
	format = Documenter.HTML(),
	modules = [LatLon]
	)

	# Documenter can also automatically deploy documentation to gh-pages.
	# See "Hosting Documentation" and deploydocs() in the Documenter manual
	# for more information.
	deploydocs(
		repo = "github.com/scottrsm/LatLon.jl.git"
	)
