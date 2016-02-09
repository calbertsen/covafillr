
.onLoad <- function(libname,pkgname){
    library.dynam("covafillr", pkgname, libname)
}

.onUnload <- function(libpath){
    library.dynam.unload("covafillr", libpath)
}
