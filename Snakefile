"""Snakemake workflow for extreme heat and drought occurrence

Note that the isimip data server appears to refuse greedy downloading,
so we can limit the intensity using a custom "download_streams" resource.
See https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources
for more details.

For example, run with:
    snakemake --cores 16 --resources download_streams=1 heat
    snakemake --cores 16 --resources download_streams=1 data/lange2020_hwmid-humidex_hadgem2-es_ewembi_rcp60_nosoc_co2_leh_global_annual_2006_2099_2030_exposure.tif

"""
from pathlib import Path
import pandas


rule all:
    input:
        "lange2020_expected_occurrence.zip",

#
# Run annual occurrence/exposure over all models/scenarios
#
rule heat:
    input:
        expand(
            "data/lange2020_hwmid-humidex_{GCM}_ewembi_{RCP}_nosoc_co2_leh_global_annual_2006_2099_{EPOCH}_{METRIC}.tif",
            GCM=["gfdl-esm2m","hadgem2-es","ipsl-cm5a-lr", "miroc5"],
            RCP=["rcp26","rcp60"],
            EPOCH=["2030","2050","2080"],
            METRIC=["occurrence", "exposure"]),
        expand(
            "data/lange2020_hwmid-humidex_{GCM}_ewembi_historical_nosoc_co2_leh_global_annual_1861_2005_{EPOCH}_{METRIC}.tif",
            GCM=["gfdl-esm2m","hadgem2-es","ipsl-cm5a-lr", "miroc5"],
            EPOCH=["baseline"],
            METRIC=["occurrence", "exposure"])

rule drought:
    input:
        # Most models have all GCMs, 2005soc in future
        expand(
            "data/lange2020_{MODEL}_{GCM}_ewembi_{RCP}_2005soc_co2_led_global_annual_2006_2099_{EPOCH}_{METRIC}.tif",
            MODEL=["clm45", "h08", "lpjml", "pcr-globwb", "watergap2"],
            GCM=["gfdl-esm2m","hadgem2-es","ipsl-cm5a-lr", "miroc5"],
            RCP=["rcp26","rcp60"],
            EPOCH=["2030","2050","2080"],
            METRIC=["occurrence", "exposure"]),

        # "clm45", "mpi-hm" use "2005soc" for baseline
        expand(
            "data/lange2020_{MODEL}_{GCM}_ewembi_historical_2005soc_co2_led_global_annual_1861_2005_{EPOCH}_{METRIC}.tif",
            MODEL=["clm45"],
            GCM=["gfdl-esm2m","hadgem2-es","ipsl-cm5a-lr", "miroc5"],
            EPOCH=["baseline"],
            METRIC=["occurrence", "exposure"]),

        # "h08", "lpjml", "pcr-globwb", "watergap2" use "histsoc" for baseline
        expand(
            "data/lange2020_{MODEL}_{GCM}_ewembi_historical_histsoc_co2_led_global_annual_1861_2005_{EPOCH}_{METRIC}.tif",
            MODEL=["h08", "lpjml", "pcr-globwb", "watergap2"],
            GCM=["gfdl-esm2m","hadgem2-es","ipsl-cm5a-lr", "miroc5"],
            EPOCH=["baseline"],
            METRIC=["occurrence", "exposure"]),

        # MPI-HM has no HadGEM
        expand(
            "data/lange2020_{MODEL}_{GCM}_ewembi_{RCP}_2005soc_co2_led_global_annual_2006_2099_{EPOCH}_{METRIC}.tif",
            MODEL=["mpi-hm"],
            GCM=["gfdl-esm2m","ipsl-cm5a-lr", "miroc5"],
            RCP=["rcp26","rcp60"],
            EPOCH=["2030","2050","2080"],
            METRIC=["occurrence", "exposure"]),
        expand(
            "data/lange2020_{MODEL}_{GCM}_ewembi_historical_histsoc_co2_led_global_annual_1861_2005_{EPOCH}_{METRIC}.tif",
            MODEL=["mpi-hm"],
            GCM=["gfdl-esm2m","ipsl-cm5a-lr", "miroc5"],
            EPOCH=["baseline"],
            METRIC=["occurrence", "exposure"]),

        # Jules-W1 and Orchidee use "nosoc"
        expand(
            "data/lange2020_{MODEL}_{GCM}_ewembi_{RCP}_nosoc_co2_led_global_annual_2006_2099_{EPOCH}_{METRIC}.tif",
            MODEL=["jules-w1", "orchidee"],
            GCM=["gfdl-esm2m","hadgem2-es","ipsl-cm5a-lr", "miroc5"],
            RCP=["rcp26","rcp60"],
            EPOCH=["2030","2050","2080"],
            METRIC=["occurrence", "exposure"]),
        expand(
            "data/lange2020_{MODEL}_{GCM}_ewembi_historical_nosoc_co2_led_global_annual_1861_2005_{EPOCH}_{METRIC}.tif",
            MODEL=["jules-w1", "orchidee"],
            GCM=["gfdl-esm2m","hadgem2-es","ipsl-cm5a-lr", "miroc5"],
            EPOCH=["baseline"],
            METRIC=["occurrence", "exposure"])

#
# Download ISIMIP extreme heat and drought timeseries/lat/lon data
#
def url_epoch(wildcards):
    if wildcards.Y == "2006":
        return "future"
    else:
        return "historical"

rule download:
    output: "incoming_data/lange2020_{MODEL}_{GCM}_ewembi_{RCP}_{SOC}_{SEN}_{VAR}_global_annual_{Y}_{Z}.nc4"
    resources: download_streams=1
    params:
        url_epoch=url_epoch
    shell:
        """
        sleep 2 && \
        wget -w 2 -nc -P ./incoming_data \
        https://files.isimip.org/ISIMIP2b/DerivedOutputData/Lange2020/{wildcards.MODEL}/{wildcards.GCM}/{params.url_epoch}/lange2020_{wildcards.MODEL}_{wildcards.GCM}_ewembi_{wildcards.RCP}_{wildcards.SOC}_{wildcards.SEN}_{wildcards.VAR}_global_annual_{wildcards.Y}_{wildcards.Z}.nc4
        """

#
# Download JRC GHSL population data
#
rule download_population_all:
    input:
        expand(
            "incoming_data/GHS_POP_E{EPOCH}_GLOBE_R2023A_{RESOLUTION}_V1_0.tif",
            EPOCH=["2020"],  # Available in: 2030, 2025, 2020, 2015, 2010, 2005, 2000, 1995, 1990, 1985, 1980, 1975
            RESOLUTION=["4326_30ss"],  # Available in: 4326_3ss, 4326_30ss, 54009_100, 54009_1000
        ),

rule download_population:
    output: "incoming_data/GHS_POP_E{EPOCH}_GLOBE_R2023A_{RESOLUTION}_V1_0.zip"
    shell:
        """
        wget -nc -P ./incoming_data \
        https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E{wildcards.EPOCH}_GLOBE_R2023A_{wildcards.RESOLUTION}/V1-0/GHS_POP_E{wildcards.EPOCH}_GLOBE_R2023A_{wildcards.RESOLUTION}_V1_0.zip
        """

rule extract_population:
    input: "incoming_data/GHS_POP_E{EPOCH}_GLOBE_R2023A_{RESOLUTION}_V1_0.zip"
    output: "incoming_data/GHS_POP_E{EPOCH}_GLOBE_R2023A_{RESOLUTION}_V1_0.tif"
    shell:
        """
        unzip {input} -d ./incoming_data
        """

#
# Find expected annual occurrence
#
EPOCH_BINS = {
    "baseline": {
        "bin_start": 1966,
        "bin_end": 2005,
        "file_year_start": 1861,
        "file_year_end": 2005,
    },
    "2030": {
        "bin_start": 2010,
        "bin_end": 2049,
        "file_year_start": 2006,
        "file_year_end": 2099,
    },
    "2050": {
        "bin_start": 2030,
        "bin_end": 2069,
        "file_year_start": 2006,
        "file_year_end": 2099,
    },
    "2080": {
        "bin_start": 2060,
        "bin_end": 2099,
        "file_year_start": 2006,
        "file_year_end": 2099,
    },
}

# e.g. heat
# data/lange2020_hwmid-humidex_gfdl-esm2m_ewembi_historical_nosoc_co2_leh_global_annual_1861_2005.nc4
# data/lange2020_hwmid-humidex_hadgem2-es_ewembi_rcp60_nosoc_co2_leh_global_annual_2006_2099_2030_occurrence.tif
rule average_nc_to_geotiff:
    input: "incoming_data/lange2020_{MODEL}_{GCM}_ewembi_{RCP}_{SOC}_{SEN}_{VAR}_global_annual_{Y}_{Z}.nc4"
    output: "data/lange2020_{MODEL}_{GCM}_ewembi_{RCP}_{SOC}_{SEN}_{VAR}_global_annual_{Y}_{Z}_{EPOCH}_occurrence.tif"
    run:
        import rioxarray
        import xarray
        ds = xarray.open_dataset(str(input), engine='netcdf4', decode_times=False)
        ds['time'] = ds.time + 1661  # correct for time in years since 1661

        epoch = EPOCH_BINS[wildcards.EPOCH]
        time_range = range(epoch["bin_start"], epoch["bin_end"] + 1)

        # find mean over range
        epoch_mean = ds.sel(time=time_range).mean(dim="time")

        # write out to TIFF
        epoch_mean.rio.write_crs("epsg:4326", inplace=True)
        epoch_mean[wildcards.VAR].rio.to_raster(str(output), compress='lzw')


rule resample:
    input:
        src_fname="incoming_data/GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif",
    output:
        dst_fname="data/GHS_POP_E2020_GLOBE_R2023A_4326_0.5deg_V1_0.tif",
    shell:
        """
        gdalwarp \
            -te -180 -90 180 90 \
            -tr 0.5 0.5 \
            -r sum \
            {input.src_fname} \
            {output.dst_fname}
        """

# e.g. heat
# data/lange2020_hwmid-humidex_hadgem2-es_ewembi_rcp60_nosoc_co2_leh_global_annual_2006_2099_2030_exposure.tif
rule population_exposure:
    input:
        population="data/GHS_POP_E2020_GLOBE_R2023A_4326_0.5deg_V1_0.tif",
        occurrence="data/lange2020_{MODEL}_{GCM}_ewembi_{RCP}_{SOC}_{SEN}_{VAR}_global_annual_{Y}_{Z}_{EPOCH}_occurrence.tif",
    output: "data/lange2020_{MODEL}_{GCM}_ewembi_{RCP}_{SOC}_{SEN}_{VAR}_global_annual_{Y}_{Z}_{EPOCH}_exposure.tif"
    shell:
        """
        gdal_calc.py \
            --calc="A*B" \
            --outfile={output} \
            -A {input.population} \
            -B {input.occurrence} \
            --co="COMPRESS=LZW" \
            --overwrite
        """

rule archive:
    input:
        rules.heat.input,
        rules.drought.input,
    output:
        archive="lange2020_expected_occurrence.zip",
    shell:
        """
        zip -r {output.archive} data
        """
