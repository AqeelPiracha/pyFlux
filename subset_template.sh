# downloading the weekly dataset of ARMOR3D from cmems servers:
# https://data.marine.copernicus.eu/product/MULTIOBS_GLO_PHY_TSUV_3D_MYNRT_015_012/services
# VARIABLE
dataset_id="$1"
variable="$2"
start_time="$3"
end_time="$4"
output_dir="$5"
echo "$start_time"
copernicusmarine subset \
	--dataset-id "$dataset_id" \
	--variable "$variable" \
	--start-datetime "$start_time" \
	--end-datetime "$end_time"  \
	--output-directory "$output_dir" \
	--force-download
