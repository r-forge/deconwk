## no S4 methodology here; speedup :
.noGenerics <- TRUE

.onLoad <- function(lib, pkg){
  version <- packageDescription("DeconWK", fields="Version")
  hello <- paste("This is DeconWK ",version,".  For an overview type 'help(\"DeconWK-package\")'.",sep="")
  packageStartupMessage(hello)
 }

.onUnload <- function(libpath)
    library.dynam.unload("DeconWK",  libpath)
