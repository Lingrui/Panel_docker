docker run --rm -it -v /usr/bin/docker:/usr/bin/docker -v /var/run/docker.sock:/var/run/docker.sock -v /usr/lib/x86_64-linux-gnu/libapparmor.so.1:/usr/lib/x86_64-linux-gnu/libapparmor.so.1 -v /nas/rawdata:/rawdata -v /nas/projects/panel/report:/report -v /nas/projects/panel/log:/log -v /nas/projects/panel/tmp:/tmp -v /ref:/ref -v /nas/common/gtf:/gtf -v /nas/common/bed:/bed \
-e RAWDATA_PATH="/nas/rawdata/" -e REPORT_PATH="/nas/projects/panel/report" -e LOG_PATH="/nas/projects/panel/log" -e TMP_PATH="/nas/projects/panel/tmp" -e REF_PATH="/ref" \
-e GTF_PATH="/nas/common/gtf" -e BED_PATH="/nas/common/bed" -e ENDPOINT=PE -e BARCODE="123" -e GENOME="/ref/hg19.fa" -e TRANSCRIPT="/gtf/hg19_gene_exon.gtf" \
-e BED="/bed/SeqCap_EZ_Exome_v3_capture.bed" -e DATASET="/nas/rawdata/rawdata_1/panel/IonXpress_012.bam"  -e \
SETNAME="control" -e SPECIFIC="no" -e THREAD="8" 192.168.6.28:5000/weixuan/runpanel:v1 bash 
