#!/bin/bash
main_dir=$1
cyto_dv_ml_dir=$2
sample=$3
size_k=$4

for svtype in manta delly cnvnator breakdancer       
    do     
       echo ${svtype}              
       svtyper-sso --core 8 --max_reads 100000 -i ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.${size_k}k -B ${main_dir}/out/${sample}/${sample}.bam > ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.${size_k}k      
#       python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_size.py ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf 10000 down > ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.${size_k}k
       python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.${size_k}k  
    done
               
# SV svtyper run    
echo "# run sytyper for all nontrs SV" && date
svtyper-sso --core 8 --max_reads 100000 -i ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all -B ${main_dir}/out/${sample}/${sample}.bam > ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all.svtyped.vcf
python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all.svtyped.vcf

echo "# run sytyper for all trs SV" && date
svtyper-sso --core 8 --max_reads 100000 -i ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all -B ${main_dir}/out/${sample}/${sample}.bam > ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all.svtyped.vcf
python ${cyto_dv_ml_dir}/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all.svtyped.vcf
 rm -rf ${main_dir}/out/${sample}/${sample}.bam*
