# Missingness quantification

Performs missing data analysis on data at different sites.

docker run \
  --name  missing4ce \
  --volume your_path_here:/4ceData \
  --rm -it \
  dbmi/4ce-analysis:version-2.1.0 R
  
remotes::install_github('covidclinical/Phase2.1MissingRPackage',
                        subdir = 'FourCePhase2.1Missing',
                        upgrade = FALSE)

library(Phase2.1MissingRPackage)
runAnalysis()


