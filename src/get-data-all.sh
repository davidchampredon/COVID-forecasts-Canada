
pids=""

for i in AB BC ON QC SK MB NS NL
do
	echo "Importing data $i ..."
	Rscript get-data.R $i &
	pids="$pids $!"
done

# Wait for all processes before continuing:
wait $pids

Rscript plot-data.R

echo "--- All Data Imported ---"
