// Earth Engine Code Editor. One model/scenario/crop year per run.
// Produces a five-band climate TIFF, a QC TIFF and a provenance CSV.
// Heat from a symmetric daily sine curve is an estimate, not hourly observations.
var YEAR = 2030;
var MODEL = 'ACCESS-CM2';
var SCENARIO = 'ssp585';
var GRID_LABEL = null; // GFDL-CM4 requires 'gr1' or 'gr2'.
var MAIZE_ASSET = 'YOUR_MAIZE_ASSET_ID';
var AREA_BAND = 'b1';
var TYPE_BAND = 'b2';
var AREA_UNIT = 'source_area_unit'; // Confirm asset units before claiming hectares.
var DRIVE_FOLDER = 'global_maize_future_cmip6_climate_v2';
var NODATA = -9999;
var BANDS = ['DD1030', 'EDH30', 'Rain'];
var REQUIRED = ['tasmin', 'tasmax', 'pr'];
var region = ee.Geometry.Rectangle([-180, -90, 180, 90], null, false);
// Inset the export request inside the outermost cells to avoid boundary padding.
// Keep the full region above for area totals and retain the product transform.
var EXPORT_INSET_DEG = 0.001;
var exportRegion = ee.Geometry.Rectangle([
  -180 + EXPORT_INSET_DEG, -90 + EXPORT_INSET_DEG,
  180 - EXPORT_INSET_DEG, 90 - EXPORT_INSET_DEG
], 'EPSG:4326', false);

if (['ssp245', 'ssp585'].indexOf(SCENARIO) < 0) {
  throw new Error('SCENARIO must be ssp245 or ssp585.');
}
if (YEAR % 1 !== 0 || YEAR < 1950 || YEAR > 2099) {
  throw new Error('YEAR must allow the following June boundary (1950..2099).');
}
if (MODEL === 'GFDL-CM4' && ['gr1', 'gr2'].indexOf(GRID_LABEL) < 0) {
  throw new Error('Choose one GFDL-CM4 grid_label to prevent double counting.');
}

// The year label is the START year of the June-May crop-year window.
var sameStart = ee.Date.fromYMD(YEAR, 1, 1);
var sameEnd = sameStart.advance(1, 'year');
var crossStart = ee.Date.fromYMD(YEAR, 6, 1);
var crossEnd = crossStart.advance(1, 'year');
var source = ee.ImageCollection('NASA/GDDP-CMIP6')
  .filter(ee.Filter.eq('model', MODEL));
if (GRID_LABEL !== null) source = source.filter(ee.Filter.eq('grid_label', GRID_LABEL));
// Stitch historical and the selected SSP at 2015 for same-method baseline runs.
var cmip6 = source.filter(ee.Filter.eq('scenario', 'historical'))
  .filterDate('1950-01-01', '2015-01-01')
  .merge(source.filter(ee.Filter.eq('scenario', SCENARIO))
    .filterDate('2015-01-01', '2101-01-01'))
  .filterDate(sameStart, crossEnd).sort('system:time_start');

function inventory(collection, start, end) {
  var subset = collection.filterDate(start, end).map(function(img) {
    var missing = ee.List(REQUIRED).removeAll(img.bandNames()).size();
    return img.set('missing_required_bands', missing)
      .set('date_key', img.date().format('YYYY-MM-dd'));
  });
  return ee.Dictionary({
    start: start.format('YYYY-MM-dd'), end_exclusive: end.format('YYYY-MM-dd'),
    expected: end.difference(start, 'day'), images: subset.size(),
    unique_dates: subset.aggregate_count_distinct('date_key'),
    missing_band_images: subset.filter(ee.Filter.gt('missing_required_bands', 0)).size(),
    provider_filled_days: subset.filter(ee.Filter.eq('interpolated', 'true')).size(),
    versions: subset.aggregate_histogram('version')
  });
}

var maize = ee.Image(MAIZE_ASSET);
// Only small metadata are retrieved. A failed preflight creates no tasks.
var check = ee.Dictionary({
  same: inventory(cmip6, sameStart, sameEnd),
  cross: inventory(cmip6, crossStart, crossEnd),
  maize_bands: maize.bandNames()
}).getInfo();
['same', 'cross'].forEach(function(key) {
  var c = check[key];
  if (c.images !== c.expected || c.unique_dates !== c.expected || c.missing_band_images !== 0) {
    throw new Error('Daily inventory failed: ' + key + ' ' + JSON.stringify(c));
  }
});
[AREA_BAND, TYPE_BAND].forEach(function(b) {
  if (check.maize_bands.indexOf(b) < 0) throw new Error('Missing maize band: ' + b);
});
print('Daily inventory (provider-filled dates retained and counted)', check);

var projection = ee.Image(cmip6.first()).select('tasmax').projection();
var grid = projection.getInfo();
if (grid.crs !== 'EPSG:4326' || Math.abs(grid.transform[0] - 0.25) > 1e-8 ||
    Math.abs(grid.transform[4] + 0.25) > 1e-8 ||
    grid.transform[1] !== 0 || grid.transform[3] !== 0) {
  throw new Error('Unexpected product grid: ' + JSON.stringify(grid));
}
var area = maize.select(AREA_BAND);
var sourceProjection = area.projection();
var cropType = maize.select(TYPE_BAND);
var positive = area.gt(0);
var validType = cropType.eq(0).or(cropType.eq(1)).unmask(0);
var validArea = area.updateMask(positive.and(validType));

// reduceResolution weights are overlap / TARGET area, not overlap / SOURCE area.
// Aggregate density, including zero background, then multiply by target area.
// Units cancel internally; output retains the input area unit.
function aggregateArea(inputArea) {
  var sourcePixelArea = ee.Image.pixelArea().reproject(sourceProjection);
  var density = inputArea.unmask(0, false).divide(sourcePixelArea)
    .setDefaultProjection(sourceProjection);
  return density.reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 1024})
    .reproject(projection).multiply(ee.Image.pixelArea().reproject(projection));
}
var area0 = aggregateArea(validArea.updateMask(cropType.eq(0)));
var area1 = aggregateArea(validArea.updateMask(cropType.eq(1)));
var nativeArea = area0.add(area1).rename('maize_area');
var maizeMask = nativeArea.gt(0);
var nativeType = area1.gt(area0).rename('crop_year_type').updateMask(maizeMask);
var minority = area0.min(area1).divide(nativeArea)
  .rename('minority_crop_year_area_fraction').updateMask(maizeMask);

function dailyClimateMetrics(img) {
  var lo = img.select('tasmin').subtract(273.15);
  var hi = img.select('tasmax').subtract(273.15);
  var rain = img.select('pr');
  var valid = img.select(REQUIRED).mask().reduce(ee.Reducer.min())
    .and(hi.gte(lo)).and(rain.gte(0)).and(maizeMask);
  var mean = lo.add(hi).divide(2);
  var amplitude = hi.subtract(lo).divide(2);
  // Numerical 24-point quadrature of the same sine curve as the original.
  var heat = ee.ImageCollection.fromImages(ee.List.sequence(0, 23).map(function(h) {
    var cosine = ee.Number(h).subtract(14).multiply(2 * Math.PI / 24).cos();
    var t = mean.add(amplitude.multiply(cosine));
    return t.min(30).subtract(10).max(0).divide(24).rename('DD1030')
      .addBands(t.subtract(30).max(0).rename('EDH30'));
  })).sum().setDefaultProjection(projection);
  return heat.addBands(rain.multiply(86400).rename('Rain'))
    .updateMask(valid).copyProperties(img, ['system:time_start']);
}
var daily = cmip6.map(dailyClimateMetrics);
function accumulate(start, end, type) {
  var subset = daily.filterDate(start, end);
  var count = subset.select('DD1030').count().unmask(0)
    .setDefaultProjection(projection);
  var typeMask = maizeMask.and(nativeType.eq(type));
  var complete = count.eq(end.difference(start, 'day')).and(typeMask);
  return {
    climate: subset.sum().setDefaultProjection(projection).updateMask(complete),
    count: count.rename('valid_days').updateMask(typeMask),
    complete: complete.rename('complete_year').updateMask(typeMask)
  };
}
var same = accumulate(sameStart, sameEnd, 0);
var cross = accumulate(crossStart, crossEnd, 1);
function combine(a, b) {
  return ee.ImageCollection.fromImages([a, b]).mosaic().setDefaultProjection(projection);
}
var climate = combine(same.climate, cross.climate);
var complete = combine(same.complete, cross.complete);
var output = climate.select(BANDS)
  .addBands(nativeArea.updateMask(maizeMask)).addBands(nativeType)
  .updateMask(complete.eq(1)).toFloat();
var qc = nativeArea.updateMask(maizeMask).addBands(nativeType).addBands(minority)
  .addBands(combine(same.count, cross.count)).addBands(complete).toFloat();

// Global source/native area checks, independent of climate missingness.
function areaTotal(img, proj) {
  return ee.Image(img).rename('area').reduceRegion({
    reducer: ee.Reducer.sum(), geometry: region,
    crs: proj, maxPixels: 1e8, tileScale: 4
  }).get('area');
}
var sourceTotal = ee.Number(areaTotal(validArea.unmask(0, false), sourceProjection));
var targetTotal = ee.Number(areaTotal(nativeArea, projection));
var invalidTypeArea = area.updateMask(positive.and(validType.not())).unmask(0, false);
var manifest = ee.Feature(null, {
  script_version: '2.1', dataset: 'NASA/GDDP-CMIP6', model: MODEL,
  export_inset_deg: EXPORT_INSET_DEG, expected_width: 1440, expected_height: 720,
  scenario_after_2014: SCENARIO, grid_label: GRID_LABEL || 'not_applicable',
  year: YEAR, maize_asset: MAIZE_ASSET, maize_mask_year: 2020, area_unit: AREA_UNIT,
  crs: grid.crs, transform: JSON.stringify(grid.transform),
  same_start: check.same.start, same_end_exclusive: check.same.end_exclusive,
  cross_start: check.cross.start, cross_end_exclusive: check.cross.end_exclusive,
  expected_days_same: check.same.expected, expected_days_cross: check.cross.expected,
  provider_filled_days_same: check.same.provider_filled_days,
  provider_filled_days_cross: check.cross.provider_filled_days,
  versions_same: JSON.stringify(check.same.versions), versions_cross: JSON.stringify(check.cross.versions),
  heat_method: '24_point_symmetric_sine_Tmin_Tmax_NOT_VALIDATED_AGAINST_ERA5',
  DD1030_unit: 'degC day', EDH30_unit: 'degC hour', Rain_unit: 'mm',
  crop_year_assignment: 'area_dominant_type_tie_to_0_start_year_label',
  source_valid_maize_area: sourceTotal, native_valid_maize_area: targetTotal,
  area_relative_difference: targetTotal.subtract(sourceTotal).divide(sourceTotal.max(1e-12)),
  excluded_invalid_crop_type_area: areaTotal(invalidTypeArea, sourceProjection),
  complete_climate_maize_area: areaTotal(nativeArea.updateMask(complete.eq(1)), projection),
  minority_crop_year_area: areaTotal(area0.min(area1), projection),
  domain: 'global_maize_requires_local_Figure2A_12_region_filter_before_RF'
});
var stem = 'MaizeFutureClimate_v2_' + MODEL + '_' + SCENARIO +
  (GRID_LABEL ? '_' + GRID_LABEL : '') + '_' + YEAR;
function exportImage(img, suffix) {
  Export.image.toDrive({
    image: img.unmask(NODATA, false), description: stem + suffix,
    folder: DRIVE_FOLDER, fileNamePrefix: stem + suffix, region: exportRegion,
    crs: grid.crs, crsTransform: grid.transform, maxPixels: 1e8,
    fileFormat: 'GeoTIFF', formatOptions: {cloudOptimized: true, noData: NODATA}
  });
}
exportImage(output, '');
exportImage(qc, '_QC');
Export.table.toDrive({collection: ee.FeatureCollection([manifest]),
  description: stem + '_manifest', folder: DRIVE_FOLDER,
  fileNamePrefix: stem + '_manifest', fileFormat: 'CSV'});
print('Product grid', grid);
print('Expected TIFF dimensions: 1440 columns x 720 rows; verify after download.');
print('Main bands', output.bandNames());
print('QC bands', qc.bandNames());
print('Run all three tasks for the test year. Review QC/manifest before batch export.');
Map.setCenter(0, 15, 2);
Map.addLayer(qc.select('minority_crop_year_area_fraction'),
  {min: 0, max: 0.5}, 'Minority crop-year fraction', false);
Map.addLayer(output.select('EDH30'), {min: 0, max: 1500}, 'Estimated EDH30', false);
