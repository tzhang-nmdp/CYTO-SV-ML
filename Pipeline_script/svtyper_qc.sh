#!/bin/bash
main_dir=$1
cyto_sv_ml_dir=$2
sample=$3
sv_caller=$4
size_k=$5

echo ${sv_caller} | sed "s%@%\n%g" > sv_caller_vector.tmp
sc_ln=$(wc -l sv_caller_vector.tmp | awk '{print $1}')
echo "# run sytyper for all sv callers" && date
for i in $(seq 1 $sc_ln)
   do
       svtype=$(awk -v a="$i" '(FNR==a){print $1}' sv_caller_vector.tmp)  
       echo ${svtype}              
       svtyper-sso --core 8 --max_reads 100000 -i ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.${size_k}k -B ${main_dir}/in/${sample}.bam > ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.${size_k}k      
       python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.${size_k}k  
    done
sudo rm -rf sv_caller_vector.tmp 

# SV svtyper run    
echo "# run sytyper for all nontrs SV" && date
svtyper-sso --core 8 --max_reads 100000 -i ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all -B ${main_dir}/in/${sample}.bam > ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all.svtyped.vcf
python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all.svtyped.vcf

echo "# run sytyper for all trs SV" && date
svtyper-sso --core 8 --max_reads 100000 -i ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all -B ${main_dir}/in/${sample}.bam > ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all.svtyped.vcf
python ${cyto_sv_ml_dir}/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all.svtyped.vcf
 rm -rf ${main_dir}/out/${sample}/${sample}.bam*
