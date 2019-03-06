
.onLoad <- function(libname = find.package("PedGFLMM"), pkgname = "PedGFLMM") {


    # CRAN Note avoidance
    if(getRversion() >= "2.15.1")
       utils::globalVariables(
         # global tables
           c("ENV")
                              )

    invisible()
}

