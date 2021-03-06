
#' plotFigure
#'
#' Recreates any figure from scratch.
#'
#' Include description.
#'
#' @name plotFigure
#' @rdname plotFigure
#' @aliases plotFigure
#' @param figure Character; Can be "go_analysis" or "RNAseq"
#' @param ... additional arguments to pass on.
#' @return spCounts object.
#' @author Jason T. Serviss
#' @examples
#'
#' \dontrun{plotFigure("go_analysis")}
#'
NULL

#' @rdname plotFigure
#' @export
#' @importFrom liftr lift render_docker
#' @importFrom tools file_path_sans_ext

plotFigure <- function(figure, ...) {
  if(figure %in% c("go_analysis", "RNAseq")) {
    path <- if_else(
      figure == "go_analysis",
      "go_analysis/go_analysis.Rmd",
      "rna_seq_analysis/rna_seq_analysis.Rmd"
    )
    runDockerAndView(path)
  } else {
    stop("You have entered an invalid figure. Please enter \"go_analysis\" or \"RNAseq\".")
  }
}

#runs liftr render_docker and opens the html output
runDockerAndView <- function(path, ...) {
  sans_ext <- file_path_sans_ext
  
  #move all files to tmp folder
  rmdPath <- system.file(path, package = "acidAdaptedRNAseq")
  tmpPath <- moveToTmp(rmdPath, ...)
  
  #render rmd in docker environment
  render_docker(file.path(tmpPath, basename(rmdPath)), cache = FALSE)
  print(tmpPath)
  
  #open resulting html
  htmlPath <- file.path(tmpPath, paste0(sans_ext(basename(rmdPath)), '.html'))
  browseURL(paste0('file://', htmlPath))
}

#runs lift and copies the .rmd file and Dockerfile to a tmp directory (due to
# the fact that liftr wants everything in the same directory)
moveToTmp <- function(rmdPath, ...){
  tmpPath <- tempdir()
  
  #copy rmd
  sysCmd1 <- paste("cp", rmdPath, tmpPath, sep = " ")
  system(sysCmd1)
  
  #lift
  lift(
    input = system.file('docker/dummy.Rmd', package = "acidAdaptedRNAseq"),
    use_config = TRUE,
    config_file = 'liftr.yml',
    output_dir = system.file('docker', package = "acidAdaptedRNAseq")
  )
  
  #copy docker
  dockerPath <- system.file('docker/Dockerfile', package = "acidAdaptedRNAseq")
  sysCmd2 <- paste("cp", dockerPath, tmpPath)
  system(sysCmd2)
  
  return(tmpPath)
}
