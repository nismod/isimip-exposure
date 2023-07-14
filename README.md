# ISIMIP Extreme Heat and Drought Occurrence and Exposure

This repository contains a data processing pipeline to produce maps of the
annual probability of drought (soil moisture below a baseline threshold) or
extreme heat (temperature and humidity-based indicators over a threshold)
events.

- Spatial resolution: 0.5° grid
- Model variables:
  - 8 hydrological models
  - 4 GCMs
  - baseline, RCP 2.6 and RCP 6.0 emission scenarios
  - current (baseline) and future maps for 2030, 2050 and 2080

Lange et al (2020) provide a timeseries of extreme events, which has been
processed into an annual probability of occurrence by the authors of this
repository.

Event definitions are given in Lange et al, Table 1:

- Land area is exposed to drought if monthly soil moisture falls below the 2.5th
  percentile of the preindustrial baseline distribution for at least seven
  consecutive months.
- Land area is exposed to extreme heat if both a relative indicator based on
  temperature (Russo et al 2015, 2017) and an absolute indicator based on
  temperature and relative humidity (Masterton & Richardson, 1979) exceed their
  respective threshold values.

This is a draft dataset, used for visualisation in
https://global.infrastructureresilience.org/ but not otherwise reviewed or
published.

If you use this, please cite:

> Lange, S., Volkholz, J., Geiger, T., Zhao, F., Vega, I., Veldkamp, T., et al.
> (2020). Projecting exposure to extreme climate impact events across six event
> categories and three spatial scales. Earth's Future, 8, e2020EF001616. DOI
> 10.1029/2020EF001616

Data citation:

> Stefan Lange, Jan Volkholz, Tobias Geiger, Fang Zhao, Iliusi Vega del Valle,
> Ted Veldkamp, Christopher Reyer, Lila Warszawski, Veronika Huber, Jonas
> Jägermeyr, Jacob Schewe, David N. Bresch, Matthias Büchner, Jinfeng Chang,
> Philippe Ciais, Marie Dury, Kerry Emanuel, Christian Folberth, Dieter Gerten,
> Simon N. Gosling, Manolis Grillakis, Naota Hanasaki, Alexandra‐Jane Henrot,
> Thomas Hickler, Yasushi Honda, Akihiko Ito, Nikolay Khabarov, Aristeidis
> Koutroulis, Wenfeng Liu, Christoph Müller, Kazuya Nishina, Sebastian Ostberg,
> Hannes Müller Schmied, Sonia I. Seneviratne, Tobias Stacke, Jörg Steinkamp,
> Wim Thiery, Yoshihide Wada, Sven Willner, Hong Yang, Minoru Yoshikawa, Chao
> Yue, Katja Frieler (2020): Land area fractions and population fractions
> exposed to extreme climate impact events derived from ISIMIP2b output data
> (v1.0). ISIMIP Repository. https://doi.org/10.48364/ISIMIP.924045

This is shared under a CC0 1.0 Universal Public Domain Dedication (CC0 1.0)

When using ISIMIP data for your research, please appropriately credit the data
providers, e.g. either by citing the DOI for the dataset, or by appropriate
acknowledgment.

The ISIMIP2b climate input data and impact model output data analyzed in this
study are available in the ISIMIP data repository at ESGF, see
https://esg.pik-potsdam.de/search/isimip/?project=ISIMIP2b&product=input and
https://esg.pik-potsdam.de/search/isimip/?project=ISIMIP2b&product=output,
respectively. More information about the GHM, GGCM, and GVM output data is
provided by Gosling et al. (2020), Arneth et al. (2020), and Reyer et al.
(2019), respectively.

## Workflow

The data processing pipeline is defined in the `Snakefile`, which uses Python,
[snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and other
dependencies which are listed in `environment.yml`.

For example, using the
[micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
python package manager (which is a smaller, faster version of the perhaps more
familiar conda or mamba):

```bash
# Install the python environment
micromamba install -f environment.yml
# Activate the python environment
micromamba activate isimip-exposure
# Run the workflow
snakemake --verbose -c32
```

## License

This code is released as open source under the MIT License, (c) 2023 Tom Russell
and contributors.

## Acknowledgments

This research received funding from the FCDO Climate Compatible Growth
Programme. The views expressed here do not necessarily reflect the UK
government's official policies.
