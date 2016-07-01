for prefix in `find *.SNP_statistics.txt | awk '{ sub(/\.SNP_statistics.txt$/, "", $0); print $0; }' | sort -u`
do 
	if [ ! -e summary_SNP_variant_statistics.xls ];then 
		sed -n '3,4p' ${prefix}.SNP_statistics.txt >>summary_SNP_variant_statistics.xls
	else sed -n '4p' ${prefix}.SNP_statistics.txt >>summary_SNP_variant_statistics.xls
	fi
	
	if [ ! -e summary_SNP_exonic_statistics.xls ];then 
		sed -n '10,11p' ${prefix}.SNP_statistics.txt >>summary_SNP_exonic_statistics.xls
	else sed -n '11p' ${prefix}.SNP_statistics.txt >>summary_SNP_exonic_statistics.xls
	fi
	
	if [ ! -e SNV.tmp1.txt ];then
		sed -n '16,17p' ${prefix}.SNP_statistics.txt >>SNV.tmp1.txt
	else sed -n '17p' ${prefix}.SNP_statistics.txt >>SNV.tmp1.txt
	fi
done

for prefix in `find *.InDel_statistics.txt | awk '{ sub(/\.InDel_statistics.txt$/, "", $0); print $0; }' | sort -u`
do 
	if [ ! -e summary_InDel_variant_statistics.xls ];then 
		sed -n '3,4p' ${prefix}.InDel_statistics.txt >>summary_InDel_variant_statistics.xls
	else sed -n '4p' ${prefix}.InDel_statistics.txt >>summary_InDel_variant_statistics.xls
	fi
	
	if [ ! -e summary_InDel_exonic_statistics.xls ];then 
		sed -n '10,11p' ${prefix}.InDel_statistics.txt >>summary_InDel_exonic_statistics.xls
	else sed -n '11p' ${prefix}.InDel_statistics.txt >>summary_InDel_exonic_statistics.xls
	fi
	
	if [ ! -e InDel.tmp1.txt ];then
		sed -n '18,19p' ${prefix}.InDel_statistics.txt >>InDel.tmp1.txt
	else sed -n '19p' ${prefix}.InDel_statistics.txt >>InDel.tmp1.txt
	fi
done

for prefix in `find *_QC.txt | awk '{ sub(/\_QC.txt$/, "", $0); print $0; }' | sort -u`
do 
	if [ ! -e summary_QC.xls ];then
		sed -n '1,2p' ${prefix}_QC.txt >>summary_QC.xls
	else sed -n '2p' ${prefix}_QC.txt >>summary_QC.xls
	fi
done

for prefix in `find *_SNP_tstv.txt | awk '{ sub(/\_SNP_tstv.txt$/, "", $0); print $0; }' | sort -u`
do 
	if [ ! -e summary_tstv.xls ];then 
		sed -n '11,12p' ${prefix}_SNP_tstv.txt >>summary_tstv.xls
	else sed -n '12p' ${prefix}_SNP_tstv.txt >>summary_tstv.xls
	fi 
	
	if [ ! -e SNV.tmp2.txt ];then 
		sed -n '5,6p' ${prefix}_SNP_tstv.txt >>SNV.tmp2.txt
	else sed -n '6p' ${prefix}_SNP_tstv.txt >>SNV.tmp2.txt
	fi 
done

for prefix in `find *_InDel_tstv.txt | awk '{ sub(/\_InDel_tstv.txt$/, "", $0); print $0; }' | sort -u`
do 
	if [ ! -e InDel.tmp2.txt ];then 
		sed -n '5,6p' ${prefix}_InDel_tstv.txt >>InDel.tmp2.txt
	else sed -n '6p' ${prefix}_InDel_tstv.txt >>InDel.tmp2.txt
	fi 
done
paste SNV.tmp1.txt SNV.tmp2.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7}' > summary_SNP_statistics.xls
paste InDel.tmp1.txt InDel.tmp2.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7}' > summary_InDel_statistics.xls
rm *.tmp*.txt


for prefix in `find *.information.xls | awk '{ sub(/\.information.xls$/, "", $0); print $0; }' | sort -u`
do 
	if [ ! -e summary.information.xls ];then 
		awk '{print $1"\t"$2}' ${prefix}.information.xls > summary.information.xls
	elif [ $prefix != summary ];then
		awk '{print $2}' ${prefix}.information.xls > temp.txt
		paste summary.information.xls temp.txt > temp2.txt
	rm temp.txt
	mv temp2.txt summary.information.xls
	fi
done

for prefix in `find *_table.txt | awk '{ sub(/\_table.txt$/, "", $0); print $0; }' | sort -u`
do
	file=${prefix##*/}
	head -6  ${prefix}_table.txt > ${file}.head_table.xls
done
