#!/bin/bash
# Example:
#	bash ./merge_simulated_samples.sh ../data/simulate_data

output_dir=$1
format=$2

if [[ -z "$output_dir" ]]; then
	echo "Usage: $0 <indir> [format: cyberDMR |HOME |Metilene |BSmooth |DSS]"
	exit 1
fi
result_dir="${output_dir}/result"
echo -e "CMD: $0 $output_dir $format"
test -d "${result_dir}" || mkdir "${result_dir}"

HEADER="chr\tstart\tend\tcoverage\tmethylation_level"
sample_names=$(find "${output_dir}" -mindepth 2 -maxdepth 2 -name "*.tsv" |xargs -n1 basename |sort -u)
for sample in ${sample_names}; do
	echo -e "Sample name is ${sample}"
	out_sample="${result_dir}/merged_${sample}"
	echo -e "$HEADER" > "${out_sample}"
	find "${output_dir}" -mindepth 2 -maxdepth 2 -name "$sample" -exec grep -vwh "start" {} + >> "${out_sample}"
	(head -n 1 "${out_sample}" && tail -n +2 "${out_sample}" |sort -k1,1V -k2,2n) > "${result_dir}/sorted_${sample}"
	rm "${out_sample}"
done
echo -e "All samples has finished"

all_formats=("cyberDMR" "HOME" "Metilene" "BSmooth" "MethyLasso" "DiffMethylTools" "DSS")

if [[ -n "$format" ]]; then
	formats=("$format")
else
	formats=("${all_formats[@]}")
fi

for format in "${formats[@]}"; do
	outdir="${result_dir}/formatted_${format}"
	test -d ${outdir} || mkdir ${outdir}

	for file in "${result_dir}"/*.tsv; do
		base=$(basename "$file")
		outfile="${outdir}/${base}"
		outfile_noheader="${outdir}/noh_${base}"

		if [[ "$format" == "cyberDMR" ]]; then
			#echo -e "Doing cyberDMR"
			echo -e "chr\tpos\tmeth_level\tcoverage" > "$outfile"
			awk 'NR>1 {print $1"\t"$2"\t"$5"\t"$4}' "$file"  >> "$outfile"
			grep -vwh "meth_level" "$outfile" > "${outfile_noheader}"
            rm "$outfile"
		elif [[ "$format" == "HOME" ]]; then
			#echo -e "Doing HOME"
			echo -e "chr\tpos\tstrand\tcontext\tmeth_count\tcoverage" > "$outfile"
			awk 'NR>1 {chr=$1; gsub("chr","",chr); meth=int($5*$4); print chr"\t"$2"\t+\tCG\t"meth"\t"$4}' "$file" >> ${outfile}
			grep -vwh "meth_count" "$outfile" > "${outfile_noheader}"
            rm "$outfile"
		elif [[ "$format" == "Metilene" ]]; then
			#echo -e "Doing Metilene"
			echo -e "chr\tstart\tend\tmethylation_level" > "$outfile"
			awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$5}' "$file"  >> "$outfile"
			grep -vwh "methylation_level" "$outfile" > "${outfile_noheader}"
            rm "$outfile"
		elif [[ "$format" == "BSmooth" ]]; then
			#echo -e "Doing BSmooth"
			echo -e "assembly\tposition\tstrand\tclass\tmc\th" > "$outfile"
			awk 'NR>1 {chr=$1; gsub("chr","",chr); meth=int($5*$4); print chr"\t"$2"\t+\tCG\t"meth"\t"$4}' "$file" >> "$outfile"
			grep -vwh "^ID" "$outfile" > "${outfile_noheader}"
            rm "$outfile"
		elif [[ "$format" == "DiffMethylTools" ]]; then
			#echo -e "Doing DiffMethylTools"
			echo -e "chr\tpos\tstrand\tpositive\tnegative\tcontxt\tcontxt2" > "$outfile"
			awk 'NR>1 {positive=int($4*$5); negative=$4-positive; print $1"\t"$2"\t+\t"positive"\t"negative"\tCG\tCGC"}' "$file" >> ${outfile}
			grep -vwh "positive" "$outfile" > "${outfile_noheader}"
            rm "$outfile"
		elif [[ "$format" == "MethyLasso" ]]; then
			#echo -e "Doing DiffMethylTools"
            cp "$file" "${outfile}"
			grep -vwh "coverage" "$outfile" > "${outfile_noheader}"
            rm "$outfile"
        elif [[ "$format" == "DSS"  ]]; then
            #echo -e "Doing DSS"
            echo -e "chr\tstart\tN\tX" > "$outfile"
            awk 'NR>1 {print $1"\t"$2"\t"$4"\t"int($4*$5)}' "$file"  >> "$outfile"
            grep -vwh "start" "$outfile" > "${outfile_noheader}"
            rm "$outfile"
		else
			echo "ERROR: Format is wrong"
			exit 2
		fi
	done
done

