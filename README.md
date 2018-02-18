# Diamond Strip Detector Analysis

Analysis software to analyse strip diamond detectors which were tested with the Strasbourg telescope


## Installation
Compile the analysis with:
```
make -j N
```
`N`: number of simultaneous jobs


## Usage
- Setup an runlist file with the runs you like to analyse according to this [template](defaultFiles/RunList.ini).
- Create run-specific settings files in a separate directory.
  The file name has to be of the form `settings.RUNNUMBER.ini` or `settings.RUNNUMBER-POSITION.ini`.
  These settings files are collected in the separate repository `StripTelescopeRunSettings`.
- Run the analysis with:
  ```
  ./diamondAnalysis -i DATA -o OUTPUT -r RUNLIST -s SETTINGS
  ```
  - `DATA`: directory of `RZ` data files
  - `OUTPUT`: output directory
  - `RUNLIST`: `ini` file to define which runs are analysed
  - `SETTINGS`: setting directory
