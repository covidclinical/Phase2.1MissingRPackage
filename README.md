# Phase2.1MissingRPackage
R code to run, validate, and submit the analysis for the Missing project.

*** We are stratifying analyses by three time periods: 
1) 1/01/2020 to 07/30/2020 ("phase_1")
2) 08/01/2020 to 09/30/2021 ("phase_2")
3) 1/01/2020 to 09/30/2021 ("phase_3")

Please run three times with each time period specified in either the runAnalysis function (for docker users) or the runAnalysis_nodocker (for non-docker users). 

To run the package:

For docker users:

```
docker run \
  --name missing4ce \
  --volume your_path_here:/4ceData \
  --rm -it \
  dbmi/4ce-analysis:version-2.4.0 R
```

While in R: 

(Change siteid to reflect your site name)

```
devtools::install_github("https://github.com/covidclinical/Phase2.1MissingRPackage", subdir="FourCePhase2.1Missing", upgrade=FALSE)
library(FourCePhase2.1Missing)

runAnalysis(dateFormat="%d-%b-%y",time = "phase_1",siteid = "penn")
submitAnalysis()
runAnalysis(dateFormat="%d-%b-%y",time = "phase_2",siteid = "penn")
submitAnalysis()
runAnalysis(dateFormat="%d-%b-%y",time = "all",siteid = "penn")
submitAnalysis()

```

For non-docker users:

```
devtools::install_github("https://github.com/covidclinical/Phase2.1MissingRPackage", subdir="FourCePhase2.1Missing", upgrade=FALSE)
library(FourCePhase2.1Missing)

data_dir = '/Your directory here /'
siteid = 'your site id'

runAnalysis_nodocker(dateFormat="%d-%b-%y",time = "phase_1")
runAnalysis(dateFormat="%d-%b-%y",time = "phase_2")
runAnalysis(dateFormat="%d-%b-%y",time = "all")
```

Please send your results files to us via slack if not using the docker image!




