

urls <- read.table("../Data_Raw/SeaIceThickness_SMOS_2015-2020/socat_search_result_url_list_20220125-210835.txt")

for (url in urls$V1) {
  download.file(url, destfile = paste0("../Data_Raw/SeaIceThickness_SMOS_2015-2020/", basename(url)))
}