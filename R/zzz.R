.First.lib <- function(libname, pkgname, where) {
	if( !require(methods) ) stop("we require methods for package: changepoint")
	
	ver <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"),"Version")
	ver <- as.character(ver)
	
	date <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"),"Date")
	date <- as.character(date)

	where <- match(paste("package:", pkgname, sep=""), search())
	.initcpt(where)
	
	cat('Successfully loaded changepoint package version ',ver,'\n')
	cat('Created on ',date,'\n')
}
