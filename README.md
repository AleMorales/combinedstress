# Stress Response Syndromes to single and combined abiotic stresses in Arabidopsis

This repository contains the scripts and data required to reproduce the paper entitled *Stress response syndromes in Arabidopsis exposed to single, simultaneous and sequential abiotic stresses*.

## Some data files missing

The raw data files from the Biosorter are not included due to Github's policies on file size and to keep the repository relatively small. However, the results of processing those files are included in the folder `Intermediate`, so the analysis can be reproduced from that point onwards. Since the processing of the raw data files is performed by the SeedSorter R package, the user may refer to that package (and associated paper) for the specifics on how the raw data from the Biosorter is processed (the same machine learning algorithms as those fitted in the original paper are used in this study).

The intermediate results from correlations are also removed for the same reason. However, this result can easily be recalculated from the existing data by running the relevant R scripts.
