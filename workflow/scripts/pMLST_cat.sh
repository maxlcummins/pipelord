#!/usr/bin/env sh
input=$(realpath $1)
output=$(realpath $2)

awk 'NR == 1 {{print "name\t" $0; next;}}{{print FILENAME "\t" $0;}}' ${input} | grep -E "pMLST profile|Sequence Type" | perl -p -e 's@pMLST profile: @@g' | awk '{{printf "%s%s",$0,NR%2?"\t":RS}}' > ${output}
#Trim the first column to generate a column with pMLST scheme
perl -p -i -e 's@^.*Inc@@g' ${output}
#Trim the first column to generate a column with pMLST scheme
perl -p -i -e 's@^.*pBSSB1@pBSSB1@g' ${output}
#Trim column two to generate a column with sample names
perl -p -i -e 's@\t.*(inc[^/]+|pbssb1[^/]+)/@\t@g' ${output}
# Trim column with sample names to clean them
perl -p -i -e 's@.out/results.txt@@g' ${output}
#Trim pMLST column to get rid of junk text
perl -p -i -e 's@Sequence Type: @@g' ${output}
#Clean square brackets from pMLST names
perl -p -i -e 's@(\]|\[)@@g' ${output}
#Fix scheme column names
perl -p -i -e 's@^@Inc@g' ${output}
perl -p -i -e 's@Inc @Inc@g' ${output}
perl -p -i -e 's@^.*family@pbssb1-family@g'  ${output}
#Change order of columns
awk -F'\t' -v OFS="\t" '{{ print $2, $1, $3}}' ${output} > tmp && mv tmp ${output}
#Insert a header
echo -e "name\tpMLST_scheme\tpMLST" | cat - ${output} > tmp && mv tmp ${output}
#Remove MLST scheme/sample combinations with no hits
grep -v 'Unknown' ${output} > tmp && mv tmp ${output}