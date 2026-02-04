# FastFrames [![build status](https://gitlab.cern.ch/atlas-amglab/fastframes/badges/main/pipeline.svg "build status")](https://gitlab.cern.ch/atlas-amglab/fastframes/commits/main) [![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14773274.svg "Zenodo DOI")](https://doi.org/10.5281/zenodo.14773274)


FastFrames is a package aimed at processing ntuples produced by [TopCPToolkit](https://topcptoolkit.docs.cern.ch/), or a similar CP algorithm based framework, into histograms or ntuples.
FastFrames rely on ROOT's [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html) to do the event loop processing.
The code allows users to define their own columns in a minimal way while using all the functionality of the processing.

## Documentation

The package documentation, including detailed instruction how to compile and run the code, can be found [here](https://cern.ch/fastframes).
Doxygen documentation for the code can be found [here](https://atlas-project-topreconstruction.web.cern.ch/fastframesdoxygen/).

## Support
Support mattermost channel is available. Please, first join the Top Analysis team: [link](https://mattermost.web.cern.ch/signup_user_complete/?id=95983da3f25882a52b0e389f0b042150&md=link&sbr=su) and then join the `Fast Frames support channel`.

## TRExFitter
FastFrames provide a convenient interface to [TRExFitter](https://gitlab.cern.ch/TRExStats/TRExFitter), a software that can be used for plotting and for running binned-likelihood fits.

## Citation
A [Zenodo record](https://zenodo.org/records/14774464) of FastFrames is available. This record is automatically updated with every version. 
We recommend citing the overall identifier: `DOI 10.5281/zenodo.14773274.`