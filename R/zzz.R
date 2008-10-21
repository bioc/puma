#.First.lib <- function(lib, pkgname, where) {
#  ## load the compiled code
#  library.dynam("puma", pkgname, lib)
#
#  if(.Platform$OS.type == "windows" && require(Biobase) && interactive()
#        && .Platform$GUI ==  "Rgui"){
#        addVigs2WinMenu("puma")
#    }
#
#  cat(paste("\nWelcome to puma version", packageDescription('puma')$Version,"\n\n"))
#
#  cat("puma is free for research purposes only. For more details, type\n")
#  cat("license.puma(). Type citation('puma') for details on how to cite\n")
#  cat("puma in publications.\n")
#}

license.puma <- function(){

  license.file <- system.file("COPYING", package="puma")

  rl <- readLines(license.file)
  rl <- rl[-c(1,length(rl))]
  rl <- gsub("*", "", rl, fixed=TRUE)
  for (i in 1:length(rl)) cat(paste(rl[i],"\n"))

}
