# ms2analyte

Data processing platform for mass spectrometry data, optimized for natural products research

## Getting Started

MS2Analyte accepts centroided mass spectrometry data from multiple instrument vendors (Waters, Thermo, Bruker) in a range of formats. mzML is the standard.

All files must have filenames in the format:

samplename_R1 where the R number is the replicate number. Samples without replicates should all be R1. Experiments with replicates should be R1, R2, R3 etc.

Samples must be in a directory named 'Samples'

Blanks (if available) must be in a directory named 'Blanks'

To run an analysis:

Clone the repository, and run pyqtgraph_test16.py

```
python pyqtgraph_test16.py
```


This starts the MS2Analyte GUI where you may select and open an existing output folder (File --> Open) or run a new analysis (Analysis --> New).

To open an existing output folder, select File and then Open from the dropdown menu. Navigate to the output folder and make sure the output folder is highlighted then click 'Open'

To open the analysis GUI, select Analysis and then New from the dropdown menu.

This starts the GUI for selecting input and output directories, and providing information about instrument and experiment variables

### Prerequisites

Requires Python 3.6+

Many dependencies. More information coming soon.

### Installing

If you do not have python already, download the latest version of python here www.python.org/downloads/

Windows users may download and install the Anaconda distribution platform here: www.anaconda.com/products/individual

Open Anaconda Prompt to create a new python environment

In the prompt, type the following:

```
conda create --name ms2analyte python=3.6
y
activate ms2analyte
```

you will then need to install the dependencies

to do that, type the following into the prompt:
```
pip install PyQt5 PyQtgraph qtmodern pandas pyteomics matchms numpy pymzml
```
Then run the program:
```
python pyqtgraph_test16.py
```
if more dependencies are needed, just use the command pip install "package" in the prompt

## Running the tests

from parent directory run:

```
python -m unittest discover
```

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/liningtonlab/ms2analyte/tags).

## Authors

* **Roger Linington** - *Initial work* - [rlinington](https://github.com/rlinington)

See also the list of [contributors](https://github.com/liningtonlab/ms2analyte/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
