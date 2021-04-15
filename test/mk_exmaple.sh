#! /bin/bash


orig_vcf=/autofs/chib/oconnor_genomes/GLAD/bing_testing_hapibd/filtered_10.vcf.gz
orig_map=/autofs/chib/oconnor_genomes/GLAD/bing_testing_hapibd/plink_map/plink.chr10.GRCh38.map


bcftools view -S <(bcftools query -l $orig_vcf | head -50) \
       	$orig_vcf \
	| bcftools annotate -x FORMAT,FILTER,INFO \
	| bcftools reheader -s <(seq 0 1 49 | sed 's/^/s/') \
	| grep -v '##bcftools' | gzip -c > test/example.vcf.gz 

cp $orig_map ./test/
mv test/"$(basename $orig_map)" test/example.map

if [ -d data ]; then mkdir data; fi
wget -P data/ https://faculty.washington.edu/browning/hap-ibd.jar
java -Xmx5G -jar data/hap-ibd.jar gt=test/example.vcf map=test/example.map out=hapibd


zcat hapibd.ibd.gz | cut -f 1,3,6,7,8 | sort -k1,1 -k2,2 -k3,3n -k4,4n \
	| awk -v OFS='\t' '{print $1 ":" $2, $3, $4, $5}' \
	> test/example.ibd

rm hapibd*.gz hapibd*.log
