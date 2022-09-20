import pandas as pd
import numpy as np
import xarray as xr
import cf_units
from datetime import datetime, timedelta
from astropy.time import Time
from dateutil.relativedelta import relativedelta

path = '/lustre/storeB/project/fou/fd/project/adc-data/tmpdata/nilu/pandora/'
file = 'groundbased_uvvis.doas.directsun.no2_nilu152_rnvs3.1.8_ny.alesund_20200313t092152z_20211002t135758z_001.nc'

unfixed = xr.open_dataset(path + file, decode_cf=False)

# fixing datetime
convert = cf_units.Unit('days since 2000-01-01 00:00:00').convert(unfixed.DATETIME.values, 'seconds since 1970-01-01 00:00:00')

new_datetime = []
for i in convert:
    new_date = datetime.fromtimestamp(i)
    new_datetime.append(new_date)
    
ds = xr.Dataset(data_vars=
                dict(altitude_instrument = (unfixed.ALTITUDE_INSTRUMENT.values),
                     angle_solar_azimuth = (['time'], unfixed.ANGLE_SOLAR_AZIMUTH.values),
                     angle_solar_zenith_astronomical = (['time'], unfixed.ANGLE_SOLAR_ZENITH_ASTRONOMICAL.values),
                     latitude = (unfixed.LATITUDE_INSTRUMENT.values),
                     longitude = (unfixed.LONGITUDE_INSTRUMENT.values),
                     no2_column_absorption_solar = (['time'], unfixed.NO2_COLUMN_ABSORPTION_SOLAR.values),
                     no2_column_absorption_solar_amf = (['time'], unfixed.NO2_COLUMN_ABSORPTION_SOLAR_AMF.values),
                     no2_column_absorption_solar_flag = (['time'], unfixed.NO2_COLUMN_ABSORPTION_SOLAR_FLAG.values),
                     no2_column_absorption_solar_uncertainty_combined_standard = (['time'], unfixed.NO2_COLUMN_ABSORPTION_SOLAR_UNCERTAINTY_COMBINED_STANDARD.values),
                     no2_column_absorption_solar_uncertainty_mixed_standard = (['time'], unfixed.NO2_COLUMN_ABSORPTION_SOLAR_UNCERTAINTY_MIXED_STANDARD.values),
                     no2_column_absorption_solar_uncertainty_random_standard = (['time'], unfixed.NO2_COLUMN_ABSORPTION_SOLAR_UNCERTAINTY_RANDOM_STANDARD.values),
                     no2_column_absorption_solar_uncertainty_systematic_standard = (['time'], unfixed.NO2_COLUMN_ABSORPTION_SOLAR_UNCERTAINTY_SYSTEMATIC_STANDARD.values),
                     temperature_effective_no2 = (['time'], unfixed.TEMPERATURE_EFFECTIVE_NO2.values),
                     temperature_effective_no2_uncertainty_combined_standard = (['time'], unfixed.TEMPERATURE_EFFECTIVE_NO2_UNCERTAINTY_COMBINED_STANDARD.values),
                     temperature_effective_no2_uncertainty_mixed_standard = (['time'], unfixed.TEMPERATURE_EFFECTIVE_NO2_UNCERTAINTY_MIXED_STANDARD.values),
                     temperature_effective_no2_uncertainty_random_standard = (['time'], unfixed.TEMPERATURE_EFFECTIVE_NO2_UNCERTAINTY_RANDOM_STANDARD.values),
                     temperature_effective_no2_uncertainty_systematic_standard = (['time'], unfixed.TEMPERATURE_EFFECTIVE_NO2_UNCERTAINTY_SYSTEMATIC_STANDARD.values)),
                coords=
                dict(time = new_datetime))

# adding variable attributes

ds['altitude_instrument'].attrs['long_name'] = 'Inst. geolocation. Altitude of the location site'
ds['altitude_instrument'].attrs['coverage_content_type'] = 'referenceInformation'
ds['altitude_instrument'].attrs['units'] = 'm'

ds['angle_solar_azimuth'].attrs['long_name'] = 'Solar azimuth for measurement center [deg], 0=north, increases clockwise'
ds['angle_solar_azimuth'].attrs['standard_name'] = 'solar_azimuth_angle'
ds['angle_solar_azimuth'].attrs['coverage_content_type'] = 'physicalMeasurement'
ds['angle_solar_azimuth'].attrs['units'] = 'degrees'

ds['angle_solar_zenith_astronomical'].attrs['long_name'] = 'Solar zenith angle for measurement center [deg]'
ds['angle_solar_zenith_astronomical'].attrs['standard_name'] = 'solar_zenith_angle'
ds['angle_solar_zenith_astronomical'].attrs['coverage_content_type'] = 'physicalMeasurement'
ds['angle_solar_zenith_astronomical'].attrs['units'] = 'degrees'

ds['latitude'].attrs['long_name'] = 'Latitude'
ds['latitude'].attrs['standard_name'] = 'latitude'
ds['latitude'].attrs['coverage_content_type'] = 'coordinate'
ds['latitude'].attrs['units'] = 'degrees_north'

ds['longitude'].attrs['long_name'] = 'Longitude'
ds['longitude'].attrs['standard_name'] = 'longitude'
ds['longitude'].attrs['coverage_content_type'] = 'coordinate'
ds['longitude'].attrs['units'] = 'degrees_east'

ds['no2_column_absorption_solar'].attrs['long_name'] = 'Nitrogen dioxide total vertical column amount [moles per square meter], -9e99=retrieval not successful'
ds['no2_column_absorption_solar'].attrs['coverage_content_type']  = 'physicalMeasurement'
ds['no2_column_absorption_solar'].attrs['units'] = 'mol m-2'

ds['no2_column_absorption_solar_amf'].attrs['long_name'] = 'Direct nitrogen dioxide air mass factor'
ds['no2_column_absorption_solar_amf'].attrs['coverage_content_type'] = 'physicalMeasurement'
ds['no2_column_absorption_solar_amf'].attrs['units'] = '1'

ds['no2_column_absorption_solar_flag'].attrs['long_name'] = 'L2 data quality flag for nitrogen dioxide, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
ds['no2_column_absorption_solar_flag'].attrs['coverage_content_type'] = 'referenceInformation'
ds['no2_column_absorption_solar_flag'].attrs['units'] = '1'

ds['no2_column_absorption_solar_uncertainty_combined_standard'].attrs['long_name'] = 'Total uncertainty of nitrogen dioxide total vertical column amount [moles per square meter]'
ds['no2_column_absorption_solar_uncertainty_combined_standard'].attrs['coverage_content_type'] = 'physicalMeasurement'
ds['no2_column_absorption_solar_uncertainty_combined_standard'].attrs['units'] = 'mol m-2'

ds['no2_column_absorption_solar_uncertainty_mixed_standard'].attrs['long_name'] = 'Structured uncertainty of nitrogen dioxide total vertical column amount [moles per square meter]'
ds['no2_column_absorption_solar_uncertainty_mixed_standard'].attrs['coverage_content_type'] = 'physicalMeasurement'
ds['no2_column_absorption_solar_uncertainty_mixed_standard'].attrs['units'] = 'mol m-2'

ds['no2_column_absorption_solar_uncertainty_random_standard'].attrs['long_name'] = 'Independent uncertainty of nitrogen dioxide total vertical column amount [moles per square meter]'
ds['no2_column_absorption_solar_uncertainty_random_standard'].attrs['coverage_content_type'] = 'physicalMeasurement'
ds['no2_column_absorption_solar_uncertainty_random_standard'].attrs['units'] = 'mol m-2'

ds['no2_column_absorption_solar_uncertainty_systematic_standard'].attrs['long_name'] = 'Common uncertainty of nitrogen dioxide total vertical column amount [moles per square meter]'
ds['no2_column_absorption_solar_uncertainty_systematic_standard'].attrs['coverage_content_type'] = 'physicalMeasurement'
ds['no2_column_absorption_solar_uncertainty_systematic_standard'].attrs['units'] = 'mol m-2'

ds['temperature_effective_no2'].attrs['long_name'] = 'Nitrogen dioxide effective temperature [K]'
ds['temperature_effective_no2'].attrs['coverage_content_type'] = 'physicalMeasurement'
ds['temperature_effective_no2'].attrs['units'] = 'degK'

ds['temperature_effective_no2_uncertainty_combined_standard'].attrs['long_name'] = 'Total uncertainty of nitrogen dioxide effective temperature [K]'
ds['temperature_effective_no2_uncertainty_combined_standard'].attrs['coverage_content_type'] = 'physicalMeasurement'
ds['temperature_effective_no2_uncertainty_combined_standard'].attrs['units'] = 'degK'

ds['temperature_effective_no2_uncertainty_mixed_standard'].attrs['long_name'] = 'Structured uncertainty of nitrogen dioxide effective temperature [K]'
ds['temperature_effective_no2_uncertainty_mixed_standard'].attrs['coverage_content_type'] = 'physicalMeasurement'
ds['temperature_effective_no2_uncertainty_mixed_standard'].attrs['units'] = 'degK'

ds['temperature_effective_no2_uncertainty_random_standard'].attrs['long_name'] = 'Independent uncertainty of nitrogen dioxide effective temperature [K]'
ds['temperature_effective_no2_uncertainty_random_standard'].attrs['coverage_content_type'] = 'physicalMeasurement'
ds['temperature_effective_no2_uncertainty_random_standard'].attrs['units'] = 'degK'

ds['temperature_effective_no2_uncertainty_systematic_standard'].attrs['long_name'] = 'Common uncertainty of nitrogen dioxide effective temperature [K]'
ds['temperature_effective_no2_uncertainty_systematic_standard'].attrs['coverage_content_type'] = 'physicalMeasurement'
ds['temperature_effective_no2_uncertainty_systematic_standard'].attrs['units'] = 'degK'

ds['time'].attrs['long_name'] = 'time'
ds['time'].attrs['standard_name'] = 'time'
ds['time'].attrs['coverage_content_type'] = 'referenceInformation'
ds['time'].attrs['units'] = 'seconds since 1970-01-01 00:00:00'
ds['time'] = ds['time'].astype(np.int32)

# Global attributes

ds.attrs['featureType'] = 'timeSeries'
ds.attrs['title'] = 'Groundbased NO2 measurements at Ny Aalesund'
ds.attrs['id'] = 'groundbased_uvvis.doas.directsun.no2_nilu152_rnvs3.1.8_ny.alesund_20200313t092152z_20211002t135758z_001'
ds.attrs['naming_authority'] = 'NILU'
ds.attrs['source'] = 'UVVIS.DOAS.DIRECTSUN.NO2_NILU152_RNVS3.1.8'
ds.attrs['summary'] = 'Remote-sensing observations performed using the Differential Optical Absorption Spectroscopy (DOAS) technique to quantify the abundance of NO2'
ds.attrs['history'] = '20220816T121440Z [geoms2nccf-0.1] groundbased_uvvis.doas.directsun.no2_nilu152_rnvs3.1.8_ny.alesund_20200313t092152z_20211002t135758z_001.h5'
ds.attrs['geospatial_lat_min'] = str(ds.latitude.values.min())
ds.attrs['geospatial_lat_max'] = str(ds.latitude.values.max())
ds.attrs['geospatial_lon_min'] = str(ds.longitude.values.min())
ds.attrs['geospatial_lon_max'] = str(ds.longitude.values.max())
ds.attrs['time_coverage_start'] = ds.time.values[0].astype('datetime64[s]').astype(datetime).strftime('%Y-%m-%d %H:%M:%S')
ds.attrs['time_coverage_end'] = ds.time.values[-1].astype('datetime64[s]').astype(datetime).strftime('%Y-%m-%d %H:%M:%S')
duration_years = str(abs(relativedelta(ds['time'].values[0].astype('datetime64[s]').astype(datetime),
                                ds['time'].values[-1].astype('datetime64[s]').astype(datetime))).years)
duration_months = str(abs(relativedelta(ds['time'].values[0].astype('datetime64[s]').astype(datetime),
                                ds['time'].values[-1].astype('datetime64[s]').astype(datetime))).months)
duration_days = str(abs(relativedelta(ds['time'].values[0].astype('datetime64[s]').astype(datetime),
                                ds['time'].values[-1].astype('datetime64[s]').astype(datetime))).days)
duration_hours = str(abs(relativedelta(ds['time'].values[0].astype('datetime64[s]').astype(datetime),
                                ds['time'].values[-1].astype('datetime64[s]').astype(datetime))).hours)
duration_minutes = str(abs(relativedelta(ds['time'].values[0].astype('datetime64[s]').astype(datetime),
                                ds['time'].values[-1].astype('datetime64[s]').astype(datetime))).minutes)
duration_seconds = str(abs(relativedelta(ds['time'].values[0].astype('datetime64[s]').astype(datetime),
                                ds['time'].values[-1].astype('datetime64[s]').astype(datetime))).seconds)
ds.attrs['time_coverage_duration'] = ('P' + duration_years + 'Y' + duration_months +
                                            'M' + duration_days + 'DT' + duration_hours + 
                                            'H' + duration_minutes + 'M' + duration_seconds + 'S')   
ds.attrs['keywords'] = 'GCMDSK:INSTRUMENTS > IN SITU/LABORATORY INSTRUMENTS > PRESSURE/HEIGHT METERS, \n' + \
                       'GCMDSK:EARTH SCIENCE > ATMOSPHERIC CHEMSITRY'
ds.attrs['keywords_vocabulary'] = 'GCMDSK:GCMD Science Keywords:https://gcmd.earthdata.nasa.gov/kms/concepts/concept_scheme/sciencekeywords,' + \
                                  'CFSTDN:CF Standard Names:https://cfconventions.org/standard-names.html'
ds.attrs['standard_name_vocabulary'] = 'CF Standard Name V79'
ds.attrs['Conventions'] = 'ACDD-1.3, CF-1.8'
ds.attrs['creator_type'] = 'Person'
ds.attrs['creator_name'] = 'Fjaeraa;Ann Mari'
ds.attrs['creator_email'] = 'ann.mari.fjaeraa@nilu.no'
ds.attrs['creator_url'] = 'nilu.no'
ds.attrs['publisher_name'] = 'Norwegian Meteorological Institute/Arcitc Data Centre'
ds.attrs['publisher_email'] = 'post@met.no'
ds.attrs['publisher_url'] = 'met.no/adc.met.no'
ds.attrs['institution'] = 'NILU'
ds.attrs['project'] = 'Arctic Data Centre'
ds.attrs['licence'] = 'https://spdx.org/licenses/CC-BY-4.0.html'
ds.attrs['references'] = 'https://www.pandonia-global-network.org/home/documents/pgn-data-use-guidelines/ \n' + \
                         'DOI: 10.48596/pgn.rnvs3p1-8.NyAlesund.P152s1'
ds.attrs['acknowledgement'] = 'The Pandonia Global Network (PGN) is a bilateral project supported with funding from NASA and ESA'

ds.to_netcdf('groundbased_no2_at_nyaalesund_20200313_20211002.nc')