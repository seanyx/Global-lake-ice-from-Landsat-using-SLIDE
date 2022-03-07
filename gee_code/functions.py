import ee

### functions to adapt to the collection1tier1 BQA band

def merge_collections_std_bandnames_collection1tier1():
  """merge landsat 5, 7, 8 collection 1 tier 1 imageCollections and standardize band names
  """
  ## standardize band names
  bn8 = ['B1', 'B2', 'B3', 'B4', 'B6', 'B7', 'BQA', 'B5', 'B10', 'B11']
  bn7 = ['B1', 'B1', 'B2', 'B3', 'B5', 'B7', 'BQA', 'B4', 'B6_VCID_1', 'B6_VCID_2']
  bn5 = ['B1', 'B1', 'B2', 'B3', 'B5', 'B7', 'BQA', 'B4', 'B6', 'B6']
  bns = ['uBlue', 'Blue', 'Green', 'Red', 'Swir1', 'Swir2', 'BQA', 'Nir', 'BT1', 'BT2']
  
  # create a merged collection from landsat 5, 7, and 8
  ls5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_TOA").select(bn5, bns)
  
  ls7 = (ee.ImageCollection("LANDSAT/LE07/C01/T1_RT_TOA")
         .filterDate('1999-01-01', '2003-05-30')
         .select(bn7, bns))
  
  ls8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_RT_TOA").select(bn8, bns)
  
  merged = ee.ImageCollection(ls5.merge(ls7).merge(ls8))
  
  return(merged)

def Unpack(bitBand, startingBit, bitWidth):
  # unpacking bit bands
  # see: https://groups.google.com/forum/#!starred/google-earth-engine-developers/iSV4LwzIW7A
  return  (ee.Image(bitBand)
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
  dem = ee.Image("users/eeProject/MERIT")##.clip(image.geometry().buffer(9000).bounds())
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


## functions to calculate ice

def CalcPre30ClimateGen():
  era5t2m = ee.ImageCollection("ECMWF/ERA5/DAILY").select(['mean_2m_air_temperature', 'total_precipitation', 'u_component_of_wind_10m', 'v_component_of_wind_10m'])
  
  def CalcPre30Climate(i):
      
    endDate = ee.Date(i.get('system:time_start'))
    startDate = endDate.advance(-30, 'day')
    mean30T2m = (era5t2m
    .filterDate(startDate, endDate)
    .mean())
    
    return ee.Image(mean30T2m)
  return CalcPre30Climate


def Assign_id_gen(i):
  def Assign_id(f):
    return f.copyProperties(i, ['LANDSAT_SCENE_ID'])
  return Assign_id

def Calc_lake_ice_gen(lakes, model, wocc):
    
  def Calc_lake_ice(i):
    snowIce = prepPredictorsTOA(i).classify(model, 'RFSnowIce').select(['RFSnowIce'])

    fmask = addFmask(i).select(['fmask'])
    cloud = fmask.eq(2).Or(fmask.eq(4)).rename(['cloud'])
    clear = cloud.Not()
    
    hillShadow = CalcHillShadowTOA(i).Not()
    
    climate = CalcPre30ClimateGen()(i)
    
    # // illustration explain sequence of masks
    # // https://docs.google.com/drawings/d/1rDYHB4aDqdAjKUQY2X2CGOmUACUzXwaF1EmQckFuyFQ/edit?usp=sharing
    i = i.addBands(fmask.eq(3).rename(['Fmask_snowIce']).updateMask(clear))
    i = i.addBands(fmask.eq(1).rename(['water']).updateMask(clear))
    i = i.addBands(fmask.eq(0).rename(['clear']).updateMask(clear))
    i = i.addBands(snowIce.rename(['SLIDE_snowIce']).updateMask(clear))
    i = i.addBands(cloud);
    
    i = (i.select(['Fmask_snowIce', 'cloud', 'water', 'clear', 'SLIDE_snowIce'])
    .updateMask(fmask.gte(0))) ## only calculate fraction when there is valid data
    
    i = i.addBands(snowIce.unmask(-1).eq(-1).rename(['missing_data'])) ## this helps with identifying cases lakes were partially outside of the landsat boundary
    
    i = i.addBands(hillShadow).addBands(climate)
    i = i.updateMask(wocc.gte(90)) ## all estimation over water surface
    
    # copyLSProperties = copyPropertiesGen(i)
    filterWithin = ee.Filter.isContained(leftField = '.geo', rightValue = i.geometry())
    lakesFil = lakes.filter(filterWithin) ##.map(copyLSProperties)

    # mergedReducer = ee.Reducer.sum().combine(ee.Reducer.mean(), None, True)
    
    pc = i.reduceRegions(
        collection = lakesFil,
        reducer = ee.Reducer.mean(),
        scale = 30,
        tileScale = 1
        ).map(Assign_id_gen(i))
    
    return pc
  return Calc_lake_ice
  
  
## lake ice model

def prepPredictorsTOA(image):
  
  # // image crs
  crs = image.select(['Nir']).projection().crs()
  
  # // scale data
  image = (image.select(['BQA'])
  .addBands(image.select(['Blue', 'Green', 'Red', 'Swir1', 'BQA', 'Nir', 'Swir2']))
  .addBands(image.select(['BT1', 'BT2']).divide(10)))
  
  # // add texture
  image = (image
  .addBands(
    image.select(['Blue', 'Green', 'Red', 'Swir1', 'Nir', 'Swir2']).multiply(100).int16()
    .glcmTexture(size = 3))
    .reproject(crs, None, 30))
    
  # // calculate the mean gradient layer
  xKernel = ee.Kernel.fixed(3, 1, [[-1, 0, 1]], 1, 0, False)
  yKernel = ee.Kernel.fixed(1, 3, [[1], [0], [-1]], 0, 1, False)
  
  gBlue = (image.select(['Blue']).convolve(xKernel)
  .addBands(image.select(['Blue']).convolve(yKernel))
  .expression('gBlue = b(0)**2 + b(1)**2').sqrt()
  .convolve(ee.Kernel.square(2, 'pixels', True))
  .reproject(crs, None, 30))
  
  image = image.addBands(gBlue)
  
  # // calculate WICI index
  WICI = (image.expression('(Swir1 + Swir2) / 2 / Nir', 
  {'Swir1': image.select('Swir1'), 'Swir2': image.select('Swir2'), 'Nir': image.select('Nir')})
  .rename('WICI'))
  
  image = image.addBands(WICI)
  
  # // add fmask
  image = addFmask(image)
  
  # // add Landsat No.
  image = (image
  .addBands(ee.Image.constant(ee.Number.parse(ee.String(image.get('LANDSAT_SCENE_ID')).slice(2, 3), 10)).select(['constant'], ['Landsat'])))
  
  # // add NDSSI
  image = image.addBands(image.normalizedDifference(['Blue', 'Nir']).rename('NDSSI'))
  
  # // add green/blue
  image = image.addBands(image.select('Green').divide(image.select('Blue')).rename('G/B'))
  
  # // add hsv
  image = image.addBands(image.select(['Red', 'Green', 'Blue']).rgbToHsv())
  
  # // add red/green
  image = image.addBands(image.select(['Red']).divide(image.select(['Green'])).rename('R/G'))
  
  return image

def toFeature(f):
  return ee.Feature(None, {'id': f})

def AddLabelGen(text):
  def AddLabel(f):
    return f.set('split', text)
  return AddLabel

def splitData(data, splitBy, seed):
  # // split data according lake ids by 70/30 
  ids = data.aggregate_array(splitBy).distinct()
  idsFC = (ee.FeatureCollection(ids.map(toFeature))
  .randomColumn('random', seed))
  trainingIds = idsFC.filterMetadata('random', 'less_than', 0.7).aggregate_array('id')
  validationIds = idsFC.filterMetadata('random', 'not_less_than', 0.7).aggregate_array('id')

  output = (data.filter(ee.Filter.inList(splitBy, trainingIds)).map(AddLabelGen('training'))
  .merge(data.filter(ee.Filter.inList(splitBy, validationIds)).map(AddLabelGen('validation'))))
  
  return output
def trainRFmodel(trainingData, ntrees, variablesPerSplit, maxNodes, inputProperties):
    
  RFclassifierTOAAll = (ee.Classifier.smileRandomForest(ntrees, variablesPerSplit, 40, 0.5, maxNodes, 2020)
  .train(
      features = trainingData, 
      classProperty = 'class_int', 
      inputProperties = inputProperties,
      subsampling = 1, 
      subsamplingSeed = 2019))
    
  return RFclassifierTOAAll

def AssignClassInt(f):
  return f.set('class_int', ee.Feature(f.get('class_int')).getNumber('class_int'))

def SLIDE(): 
  trainingToaSr = ee.FeatureCollection("users/eeProject/Lake_ice_classification/lake_ice_training_data_TOA_SR_4a18cc9d3762c9b56d7ea5fd7f631e9b")
  predictors = ee.List(['Nir', 'WICI', 'gBlue', 'Green', 'hue', 'saturation', 'value', 'G/B', 'R/G', 'NDSSI', 'Landsat'])
  
  # // assign integer class labels to classes at the same time merge certain classes
  classLabels = ['water', 'clear_ice', 'opaque_ice', 'snow']
  classLabelsTable = (ee.FeatureCollection([
                                            ee.Feature(None, {'class': 'water', 'class_int': 0, 'color': 'blue'}),
                                            ee.Feature(None, {'class': 'clear_ice', 'class_int': 1, 'color': 'orange'}),
                                            ee.Feature(None, {'class': 'opaque_ice', 'class_int': 1, 'color': 'cyan'}),
                                            ee.Feature(None, {'class': 'snow', 'class_int': 1, 'color': 'cyan'})]))
  training = (ee.Join.saveFirst(matchKey = 'class_int', outer = False).apply(trainingToaSr, classLabelsTable, ee.Filter.equals(leftField = 'class', rightField = 'class'))
  .map(AssignClassInt))
  
  trainingTOA = training.filterMetadata('source', 'equals', 'TOA').filterMetadata('WICI', 'not_equals', None)
  # // select data from which satellite to use
  # // trainingTOA = trainingTOA.filter(ee.Filter.inList('Landsat', [5, 8]));
  trainingTOA = trainingTOA.select(
    propertySelectors = predictors.cat(ee.List(['Hylak_id', 'class_int', 'fmask', 'LANDSAT_SCENE_ID'])), 
    retainGeometry = False);
  
  splitedData = splitData(trainingTOA, 'Hylak_id', 2019)
  training = splitedData.filterMetadata('split', 'equals', 'training')
  validation = splitedData.filterMetadata('split', 'equals', 'validation')
  
  model = trainRFmodel(training, 250, 2, 24, predictors)
  
  return model
