# // This script calculate lake ice percentage from landsat imageries
# // by: Xiao Yang
# // 11/04/2018

import ee
import argparse
import math
from functions import *

parser = argparse.ArgumentParser(prog = 'Download_landsat_lake_ice_data.py', description = "")

parser.add_argument('-sy', '--start_year', help = 'starting (inclusive) year of the export data', type = int, default = 2001)
parser.add_argument('-sm', '--start_month', help = 'starting (inclusive) month of the export data', type = int, default = 8)
parser.add_argument('-ey', '--end_year', help = 'end (inclusive) year of the export data', type = int, default = 2001)
parser.add_argument('-em', '--end_month', help = 'end (inclusive) month of the export data', type = int, default = 9)

parser.add_argument('-minc', '--minimal_cloud_score', help = 'minimal (inclusive) cloud score to include', type = int, default = 0)
parser.add_argument('-maxc', '--maximum_cloud_score', help = 'maximum (inclusive) cloud score to include', type = int, default = 25)

parser.add_argument('-o', '--output_dir', help = 'output directory in the Google drive', type = str, default = 'default')

parser.add_argument('--restart_no', help = 'index to start the task from', type = int, default = 0)

args = parser.parse_args()
start_year = args.start_year
start_month = args.start_month
end_year = args.end_year
end_month = args.end_month
minc = args.minimal_cloud_score
maxc = args.maximum_cloud_score
output_dir = args.output_dir
rs = args.restart_no

n = (end_year - start_year) * 12 + (end_month - start_month) + start_month + 1

ee.Initialize()

ls = (merge_collections_std_bandnames_collection1tier1()
.filterMetadata('CLOUD_COVER_LAND', 'not_less_than', minc)
.filterMetadata('CLOUD_COVER_LAND', 'not_greater_than', maxc))

lakesFil = (ee.FeatureCollection("users/eeProject/HydroLAKES_polys_v10")
.filterMetadata('Lake_area', 'greater_than', 1)
.select(['Hylak_id']))

model = SLIDE()

wocc = ee.Image("JRC/GSW1_3/GlobalSurfaceWater").select(['occurrence'])

exportedProperties = ['system:index', 'LANDSAT_SCENE_ID', 'Fmask_snowIce', 'cloud', 'water', 'clear', 'SLIDE_snowIce', 'missing_data', 'hillshadow', 'mean_2m_air_temperature', 'total_precipitation', 'u_component_of_wind_10m', 'v_component_of_wind_10m', 'Hylak_id']

if output_dir == 'default':
    output_dir = 'SLIDE_lake_ice_' + 'cs' + str(minc) + '-' + str(maxc) + '_' + str(start_year) + '-' + str(start_month) + '_' + str(end_year) + '-' + str(end_month)

for i in range(start_month + rs, n):

    yearOffset = math.floor((i - 1) / 12)
    year = int(start_year + yearOffset)
    month = int(i - yearOffset * 12)

    year_offset_next = math.floor(i / 12)
    year_next = int(start_year + year_offset_next)
    month_next = int(i + 1 - year_offset_next * 12)

    t_start = str(year) + '-' + format(month, '02d') + '-01'
    t_end = str(year_next) + '-' + format(month_next, '02d') + '-01'

    lsFil = ls.filterDate(t_start, t_end)

    monthly_lake_ice = lsFil.map(Calc_lake_ice_gen(lakesFil, model, wocc)).flatten().limit(100)
    
    # print(ee.Feature(monthly_lake_ice.first()).toDictionary().getInfo())
    fn = 'SLIDE_lake_ice' + format(minc, '02d') + '-' + format(maxc, '02d') + '_' + t_start

    task = (ee.batch.Export.table.toDrive(
            collection = monthly_lake_ice.select(propertySelectors = exportedProperties, retainGeometry = False),
            description = fn,
            folder = output_dir,
            fileNamePrefix = fn,
            fileFormat = 'CSV'))

    task.start()

    maximum_no_of_tasks(8, 120)

    print('task', i - start_month + 1, 'of', n - start_month, ', (', t_start, t_end, ')', 'has submitted.')
