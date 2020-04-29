
rm -rf fcst*.RData

pids=""

for prov in ON AB QC 
do
	for scen in BASELINE INST_RED_30 RED_50_21d
	do
		echo $prov $scen
		Rscript fcst.R $prov $scen  &
		pids="$pids $!"
	done
done

# Run analysis once all is done:
wait $pids

# Digestion and analysis:
Rscript digest.R
Rscript analyze.R > out/analyze.txt


# scen: BASELINE INST_RED_30  RED_30_21d PULSE_14d_7d_50
