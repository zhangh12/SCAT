
.onAttach <- function(libname, pkgname) {
  if (interactive()) {
    packageStartupMessage('SCAT ', packageVersion('SCAT'), '  For help type ?scat')
    packageStartupMessage('The most frequently updated version can be downloaded from https://github.com/zhangh12/SCAT')
  }
}
