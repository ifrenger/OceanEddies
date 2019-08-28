# OceanEddies
This is a  collection of algorithms developed to autonomously identify and track mesoscale ocean
eddies as closed contours of anomalies of sea surface height (SSH) or isopycnal layer depths. In principle, it can be adjusted to various applications that require an identification of closed contours of a property containing a single extremum.

# Support
The repository is forked from https://github.com/jfaghm/OceanEddies. We have used the code ([frozen version](https://zenodo.org/record/13037#.XV6C3pMzZSw), equivalent to release [v1.1](https://github.com/ifrenger/OceanEddies/releases)) in [Faghmous et al, 2015](https://www.nature.com/articles/sdata201528), where we published a [dataset]( http://dx.doi.org/10.5061/dryad.gp40h) of ocean mesoscale eddies, identified globally based on satellite sea level anomalies. I am supporting here a subsample of the code, see example application below. 

# Requirements
 + Matlab, including Matlab Mapping Toolbox

# Sample usage
The script [sample_application.m](sample_application.m) exemplifies a standard application, that is, eddies identified based on SSH anomalies. The script uses three main functions:
1. [scan_single.m](eddyscan/scan_single.m) in the directory eddyscan identifies eddies based on the input field such as SSH anomalies.
2. [tolerance_track_lnn.m](track_lnn/tolerance_track_lnn.m) tracks identified eddies over time, allowing for eddies to "disappear" for a few conecutive time steps.
3. [reformat_track_data_to_chelton.m](track_lnn/reformat_track_data_to_chelton.m) reformats data to a "matrix-style" format, similar to the one the Chelton-data set offers.

If you pull the repository the script should just run. It uses sample data provided in the data directory and produces sample plots in the output plot directory.
