# pycdata
<i>Author: Jayshil A Patel ([jayshil.patel@astro.su.se](mailto:jayshil.patel@astro.su.se))</i>

`pycdata` is a module to import datasets from various telescopes/instruments in [`pycheops`](https://github.com/pmaxted/pycheops). [`pycheops`](https://github.com/pmaxted/pycheops) is a tool that is specifically designed to model CHEOPS observations of transits, eclipses and phase curves. While being a genius tool, what it lacks is a facility to model transits (and eclipses and phase curves) from other telescopes/instruments, even the PSF photometry produced by [`PIPE`](https://github.com/alphapsa/PIPE). `pycdata` is here to fill up this very need! A good thing about `pycdata` is that it automatically put resultant data products to the [`pycheops`](https://github.com/pmaxted/pycheops) cache directory so that it can be directly readable from [`pycheops`](https://github.com/pmaxted/pycheops) command line.

Below, you will find a general instruction on its installation, usage and resultant data products.

## Installation
Installation for `pycdata` can be done using `setup.py` file in the repository, by following commands below:

```
git clone https://github.com/Jayshil/pycdata.git
cd pycdata
python setup.py install
```

There you are! You are now ready to use this package!

## PIPE Data
As already been mentioned, [`PIPE`](https://github.com/alphapsa/PIPE) produces PSF photometry for the CHEOPS targets. To make the data products from [`PIPE`](https://github.com/alphapsa/PIPE) accessible to [`pycheops`](https://github.com/pmaxted/pycheops), follow,

```
from pycdata import dpr
dpr.pipe_data(name, fileid, imagette=True)
```

`name` is the name of the fits file (data products from [`PIPE`](https://github.com/alphapsa/PIPE)), and `fileid` is a unique file key for each of the CHEOPS observations. The keyword, which is boolean, `imagette` can be set `True` if [`PIPE`](https://github.com/alphapsa/PIPE) used CHEOPS imagettes to produce PSF photometry. The columns for which the data was not available, e.g., contamination, dark, and smearing was set to zero.

## TESS Data
It is even simpler to use `pycdata` with TESS (the Transiting Exoplanets Survey Satellite) data products:

```
from pycdata import dpr
dpr.tess_data(name, pdc=True, verbose=True)
```

By providing the name of the target to the `name` keyword, and choosing whether to use PDC-SAP flux with `pdc` boolean, the TESS data products, which are readable to [`pycheops`](https://github.com/pmaxted/pycheops), can be downloaded. `pycdata` uses [`astroquery`](https://astroquery.readthedocs.io/en/latest/index.html) package to download TESS data products directly from the Mikulski Archive for Space Telescopes ([MAST](https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html)). The resultant data products contain time, pdc-sap (or sap) flux, uncertainty in flux, and centroids of the aperture. The rest of the columns (which are imporant for CHEOPS but not available for TESS, e.g., roll angle, smear etc.) are set to zero.

## Kepler/K2 Data
The procedure to download Kepler/K2 lightcurves are mostly similar to that in case of TESS data:

```
from pycdata import dpr
dpr.kepler_data(name, pdc=True, long_cadence=True, verbose=True)
```

The only difference is the `long_cadence` keyword since Kepler mission offers its data products in two cadences.

## Comments/Suggestions/Contributions
We would love to hear your comments and/or suggestions to make `pycdata` more accessible. Furthermore, if you want to make any contributions to the project, you are more than welcome --- feel free to open a pull request.