import ee

### functions to adapt to the collection1tier1 BQA band

def merge_collections_std_bandnames_collection1tier1():
    """merge landsat 5, 7, 8 collection 1 tier 1 imageCollections and standardize band names
    """
    ## standardize band names
    bn8 = ['B2', 'B3', 'B4', 'B6', 'BQA', 'B5']
    bn7 = ['B1', 'B2', 'B3', 'B5', 'BQA', 'B4']
    bn5 = ['B1', 'B2', 'B3', 'B5', 'BQA', 'B4']
    bns = ['Blue', 'Green', 'Red', 'Swir1', 'BQA', 'Nir']

    # create a merged collection from landsat 5, 7, and 8
    ls5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_TOA").select(bn5, bns)

    ls7 = (ee.ImageCollection("LANDSAT/LE07/C01/T1_RT_TOA")
           .filterDate('1999-01-01', '2003-01-01')
           .select(bn7, bns))

    ls8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_RT_TOA").select(bn8, bns)

    merged = ee.ImageCollection(ls5.merge(ls7).merge(ls8))

    return(merged)

def Unpack(bitBand, startingBit, bitWidth):
    # unpacking bit bands
    # see: https://groups.google.com/forum/#!starred/google-earth-engine-developers/iSV4LwzIW7A
    return (ee.Image(bitBand)
            .rightShift(startingBit)
            .bitwiseAnd(ee.Number(2).pow(ee.Number(bitWidth)).subtract(ee.Number(1)).int()))

def UnpackAll(bitBand, bitInfo):
    unpackedImage = ee.Image.cat([Unpack(bitBand, bitInfo[key][0], bitInfo[key][1]).rename([key]) for key in bitInfo])
    return unpackedImage

def addFmask(image):
    bitInfo = {
        'Cloud': [4, 1],
        # 'CloudConfidence': [5, 2],
        'CloudShadowConfidence': [7, 2],
        'SnowIceConfidence': [9, 2]
    }
    temp = UnpackAll(image.select(['BQA']), bitInfo)

    ## define fmask water manually
    ndvi = image.normalizedDifference(['Nir', 'Red'])
    nir = image.select(['Nir'])
    fwater = ndvi.lt(0.01).And(nir.lt(0.11)).Or(ndvi.lt(0.1).And(nir.lt(0.05)))

    fmask = (fwater.rename(['fmask'])
    .where(temp.select(['SnowIceConfidence']).eq(3), ee.Image(3))
    .where(temp.select(['CloudShadowConfidence']).eq(3), ee.Image(2))
    .where(temp.select(['Cloud']), ee.Image(4))
    ).mask(temp.select(['Cloud']).gte(0)) ## mask the fmask so that it has the same footprint as the BQA band

    return image.addBands(fmask)

def CalcHillShadowTOA(image):
    dem = ee.Image("users/eeProject/MERIT").clip(image.geometry().buffer(9000).bounds())
    SOLAR_AZIMUTH_ANGLE = ee.Number(image.get('SUN_AZIMUTH'))
    SOLAR_ZENITH_ANGLE = ee.Number(90).subtract(ee.Number(image.get('SUN_ELEVATION')))

    return(ee.Terrain.hillShadow(dem, SOLAR_AZIMUTH_ANGLE, SOLAR_ZENITH_ANGLE, 100, True)
    .reproject("EPSG:4326", None, 90).rename(['hillshadow']))

def maximum_no_of_tasks(MaxNActive, waitingPeriod):
	"""maintain a maximum number of active tasks
	"""
	import time
	import ee
	ee.Initialize()

	time.sleep(2)
	## initialize submitting jobs
	ts = list(ee.batch.Task.list())

	NActive = 0
	for task in ts:
		if ('RUNNING' in str(task) or 'READY' in str(task)):
			NActive += 1
	## wait if the number of current active tasks reach the maximum number
	## defined in MaxNActive
	while (NActive >= MaxNActive):
		time.sleep(waitingPeriod) # if reach or over maximum no. of active tasks, wait for 2min and check again
		ts = list(ee.batch.Task.list())
		NActive = 0
		for task in ts:
			if ('RUNNING' in str(task) or 'READY' in str(task)):
				NActive += 1
	return()
