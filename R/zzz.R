.onAttach <- function(libname, pkgname)
{
	f <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"),
                      c("Version", "Date"))
        packageStartupMessage('Successfully loaded changepoint package version ',
                              f[1,1],'\n', 'Created on ', f[1,2],'\n Substantial changes to the structure of the package have occured between version 0.8 and 1.0.  Please see the package NEWS for details.\n')
}
