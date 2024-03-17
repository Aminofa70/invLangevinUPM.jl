using Documenter
using invLangevinUPM

makedocs(
    sitename = "invLangevinUPM",
    format = Documenter.HTML(),
    modules = [invLangevinUPM]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
