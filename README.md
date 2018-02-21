# Diamond Strip Detector Analysis

Analysis software to analyse strip diamond detectors which were tested with the Strasbourg telescope


## Installation
Compile the analysis with:
```
make -j N
```
`N`: number of simultaneous jobs


## Usage
- Setup a runlist file with the runs you like to analyse according to this [template](defaultFiles/RunList.ini).
  The following parameters are set:

  Parameter  | Description
  ---------- | -----------
  `RunNo`    | run number
  `RunDes`   | position of the diamond sample on the readout chip (`0`, `left` or `right`)
  `Verb`     | verbosity
  `NEvents`  | number of events to be analysed
  `NStart`   | first event
  `Ped`      | calculate pedestals
  `Clus`     | clustering
  `Sel`      | apply selection
  `DoAli`    | alignment
  `AnaAli`   | analysis of alignment
  `TransAna` | transparent analysis

- Create run-specific settings files in a separate directory.
  All run parameters, such as masked channels, fiducial region, cluster thresholds and alignment options are defined in this file.
  The file name has to be of the form `settings.RUNNUMBER.ini` or `settings.RUNNUMBER-POSITION.ini`.
  A template is available [here](defaultFiles/settings.ini).
  These settings files are collected in the separate repository `StripTelescopeRunSettings`.
- Run the analysis with:
  ```
  ./diamondAnalysis -i DATA -o OUTPUT -r RUNLIST -s SETTINGS
  ```
  - `DATA`: directory of `RZ` data files
  - `OUTPUT`: output directory
  - `RUNLIST`: `ini` file to define which runs are analysed
  - `SETTINGS`: setting directory
