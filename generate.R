library(stars)
library(dplyr)
library(sf)
library(pals)

d  = dir('/Users/mjohnson/Downloads/36', full.names = TRUE)

d2 = grep('36/005|36/047|36/047|36/061|36/081|36/085', d, value = TRUE)


for(i in 1:length(d2)){
  
  f = list.files(d2[[i]], full.names = TRUE)
  out = paste0("temp/1km/36", basename(d[[i]]), '.tif')
  unlink(out)
  
  rr = arrow::open_dataset(f) %>% 
    select(X, Y) %>% 
    collect() %>% 
    st_as_sf(coords = c("X", "Y"), crs = 4326) %>% 
    st_transform(5070) %>% 
    mutate(count = 1)
  
  ss = sf::st_bbox(rr) %>% 
    stars::st_as_stars(values = 0, 
                       dx = 1000, 
                       dy = 1000)
  
  s = stars::st_rasterize(rr[, "count"], ss, options = "MERGE_ALG=ADD")
  s$count = ifelse(s$count == 0, NA, s$count)
  
  stars::write_stars(s, out, "count")
  message(i)
}

merge_rasters <- function(input_rasters,
                          output_raster = tempfile(fileext = ".tif"),
                          options = character(0)) {
  
  unlink(output_raster)
  sf::gdal_utils(
    util = "warp",
    source = as.character(input_rasters),
    destination = output_raster,
    options = options
  )
}





# rbind(c(0, 0, 0, 0),
#       cbind(1:100, t(col2rgb(rev((pals::brewer.ylorbr(100))
#       ))))) %>%
#   write.table("col.txt", row.names = FALSE, col.names = FALSE)

files = list.files('temp/1km/39', full.names = TRUE)


merge_rasters(list.files('temp/1km/', full.names = TRUE), 
              output_raster = 'temp/ohio1km.tif',
              c('-srcnodata', 'nan',
                '-dstnodata', 0,
                "-s_srs", "EPSG:5070",
                "-t_srs", "EPSG:4326",
                '-r', "near"))

system('gdal_translate -of GTiff -ot byte temp/ohio1km.tif temp/ohio1km_8.tif -scale 0 518')
system('gdaldem color-relief temp/ohio1km_8.tif col.txt temp/ohio1km_8_col.tif')
system('gdal2tiles.py --zoom=4-7 /Users/mjohnson/github/urmi-vis-ohio/ohio8_col.tif tiles1km')
