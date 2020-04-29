#
# FIT WITH ABC THEN FORECAST
#

rm -rf fcst*.RData

pids=""

for prov in ON #AB QC MB BC NL SK
do
	for scen in BASELINE ISO1 ISO2
	do
		echo $prov $scen
		Rscript fit-abc.R $prov $scen  
		pids="$pids $!"
	done
done

# Run analysis once all is done:
# wait $pids

# Digestion and analysis:
Rscript digest.R
Rscript analyze.R > out/analyze.txt



