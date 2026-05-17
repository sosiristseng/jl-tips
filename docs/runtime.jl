# # Runtime environment
import InteractiveUtils
InteractiveUtils.versioninfo()

# Direct dependencies
using Pkg
Pkg.status()

# All dependencies (manifest)
Pkg.status(mode=PKGMODE_MANIFEST)
