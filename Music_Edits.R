library(MuSiC)

# This function will get called instead of stats::cov, no need to reassign in package
cov_orig <- stats::cov
cov <- function(input, ...) {
  return(cov_orig(x = input, use = "na.or.complete", ...))
}

# This function we have to reassign in the MuSiC package in order for it to work.
# This is REALLY dirty and relies on a global environment variable called
# "sc_basis_precomputed" being defined (it isn't passed into the function), but
# also avoids me having to duplicate the entire "music_prop" function to make a
# minor change
music_basis_orig <- MuSiC::music_basis
music_basis <- function(...) {
  if (is.null(sc_basis_precomputed)) {
    print("Computing music_basis...")
    return(music_basis_orig(...))
  }
  else {
    print("Using pre-computed sc.basis...")
    print(dim(sc_basis_precomputed$Disgn.mtx))
    return(sc_basis_precomputed)
  }
}

R.utils::reassignInPackage("music_basis", pkgName="MuSiC", music_basis)

# This function we have to reassign in the MuSiC package in order for it to work.
# MuSiC will throw an error if weights are negative. This sets negative weights
# to zero.
weight.cal.ct.orig <- MuSiC::weight.cal.ct
weight.cal.ct <- function (...) {
  weight = weight.cal.ct.orig(...)
  if (any(weight < 0)) {
    print("WARNING: Weight matrix has negative weights!")
  }
  weight[weight < 0] = 0
  return(weight)
}

R.utils::reassignInPackage("weight.cal.ct", pkgName="MuSiC", weight.cal.ct)
