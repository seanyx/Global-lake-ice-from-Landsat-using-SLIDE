# // This script calculate lake ice percentage from landsat imageries
# // by: Xiao Yang
# // 11/04/2018

import ee
import argparse
import math
import numpy as np
import pandas as pd
from functions import *

parser = argparse.ArgumentParser(prog = 'Landsat_lake_ice_using_SLIDE_from_file.py', description = "")

parser.add_argument('DATE_FILE', help = 'CSV file contains yr1, yr2, mth1, mth2 columns', type = str)

parser.add_argument('-minc', '--minimal_cloud_score', help = 'minimal (inclusive) cloud score to include', type = int, default = 0)
parser.add_argument('-maxc', '--maximum_cloud_score', help = 'maximum (inclusive) cloud score to include', type = int, default = 25)

parser.add_argument('-o', '--output_dir', help = 'output directory in the Google drive', type = str, default = 'python_gee_v1_all_columns')

args = parser.parse_args()
DATE_FILE = args.DATE_FILE
minc = args.minimal_cloud_score
maxc = args.maximum_cloud_score
output_dir = args.output_dir

yrMths = pd.read_csv(DATE_FILE)
yr1 = yrMths['yr1'].values.tolist()
yr2 = yrMths['yr2'].values.tolist()
mth1 = yrMths['mth1'].values.tolist()
mth2 = yrMths['mth2'].values.tolist()

N = len(yr1)

ee.Initialize()

ls = (merge_collections_std_bandnames_collection1tier1()
.filterMetadata('CLOUD_COVER_LAND', 'not_less_than', minc)
.filterMetadata('CLOUD_COVER_LAND', 'not_greater_than', maxc))

lakesFil = (ee.FeatureCollection("users/eeProject/HydroLAKES_polys_v10")
.filterMetadata('Lake_area', 'greater_than', 1)
.select(['Hylak_id']))

model = SLIDE()

wocc = ee.Image("JRC/GSW1_3/GlobalSurfaceWater").select(['occurrence'])

exportedProperties = ['LANDSAT_SCENE_ID', 'Fmask_snowIce', 'cloud', 'water', 'clear', 'SLIDE_snowIce', 'missing_data', 'hillshadow', 'mean_2m_air_temperature', 'total_precipitation', 'u_component_of_wind_10m', 'v_component_of_wind_10m', 'Hylak_id']

for i in range(N):
    
    year = int(yr1[i])
    year_next = int(yr2[i])
    month = int(mth1[i])
    month_next = int(mth2[i])

    t_start = str(year) + '-' + format(month, '02d') + '-01'
    t_end = str(year_next) + '-' + format(month_next, '02d') + '-01'

    lsFil = ls.filterDate(t_start, t_end)

    monthly_lake_ice = lsFil.map(Calc_lake_ice_gen(lakesFil, model, wocc)).flatten()
    
    fn = 'SLIDE_lake_ice' + format(minc, '02d') + '-' + format(maxc, '02d') + '_' + t_start

    task = (ee.batch.Export.table.toDrive(
            collection = monthly_lake_ice.select(propertySelectors = exportedProperties, retainGeometry = False),
            description = fn,
            folder = output_dir,
            fileNamePrefix = fn,
            selectors = exportedProperties,
            fileFormat = 'CSV'))

    task.start()

    maximum_no_of_tasks(4, 600)
    # print(t_start, t_end)
    print(fn, 'has submitted.')
