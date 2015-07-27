.onAttach <- function(libname, pkgname)
{
	f <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"),
                      c("Version", "Date"))
        packageStartupMessage('Successfully loaded changepoint package version ',
                              f[1,1],'\n Substantial changes to the structure of the package have occurred from version 2.0. Please see the package NEWS for details.\n NOTE: Default penalty has changed to MBIC so default code may return different segmentations.')
}