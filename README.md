# EscanabaPhenology
Analysis of Escanaba Lake, WI, walleye spawning phenology.  For questions, contact zsfeiner@wisc.edu

## Scripts files
### PredictWaterTemps.R 
Uses Sparkling Lake, WI, under-ice water temperatures to build a model to predict under-ice water temps based on air temperature and freeze and thaw dates, then applies it to Escanaba data to predict Escanaba Lake under-ice water temperature data for use in recruitment modeling.

### Recruitment_Reanalysis.R S
Summarizes and cleans data and performs GAMM, HGAM, and climate windows modeling, including creation of all results, tables, and figures shown in Feiner et al. 2025.

## Data files
### EscanabaAge0PE_raw.csv
Age-0 walleye catch information from fall night electrofishing, indexes recruitment and allows for calculation of PE and PE CV
### EscanabaTemps_1956_2020.csv
Daily observed weather and water temperature data from Escanaba Lake, WI
### ModeledEscanabaWaterTemps_1956.2020.csv
Modeled Escanaba Lake water temperatures, from PredictWaterTemps.R
### NOAA_NCDC_NorthernWIWeather_1940_2020.csv
NOAA NCDC weather data for weather stations around northern WI, USA
### SparklingLakeIceDates.txt
Sparkling Lake ice-on and ice-off dates, from NTL LTER
### SparklingLakeWinterTemps.txt
Sparkling Lake water temperature data, from NTL LTER