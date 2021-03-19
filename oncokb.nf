println "***********************"
println "Yabili_Onco_Anno_Pipe"


//PARTE_1
//params.output_maf='/hpcshare/genomics/yabili_pipeline_vda/sarek_analyses/results/Annotation/10643_vs_2682/VEP/out_filtered_SOMATIC.maf'
//aggiungere un if per Strelka vs Mutect
//params.input_vcf='./results/VariantCalling/Mutect2/filt_vcf/*
//params.input_vcf="$PWD/results/VariantCalling/*_vs_*/Strelka/{StrelkaBP_*_vs_*_somatic_indels.vcf,StrelkaBP_*_vs_*_somatic_snvs.vcf}"
//StrelkaBP_*_vs_*_somatic_indels.vcf.gz
//StrelkaBP_*_vs_*_somatic_snvs.vcf.gz


params.input_vcf=""
params.tumor_id=''
params.normal_id=''
params.param=''

//fixed data
params.ref_fasta='/hpcshare/genomics/references/gatk_bundle/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta'
params.vep_path='/usr/local/bin'
params.vep_data='/hpcshare/genomics/yabili_pipeline_vda/yabili_resources_vda/vep_cache/vep_102'
params.ncbi_build='GRCh38'

//PARTE_2
IMAF="è il dato che esce da process vcf2maf"
OMAF="è il dato che esce da oncokb"
IC="data/example_clinical.txt"  
//SAMPLE_ID	ONCOTREE_CODE
//https://github.com/cBioPortal/oncotree
//https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7233083/
OC="data/example_clinical.oncokb.txt" //da aggiungere

//fixed data
TOKEN="d8bc17e4-4d95-4104-b5ee-1a2c8a94b7aa"

params.outdir="results/Annotation/OncoKB"

//Print some stuff here
println "input: $params.input_vcf"
println "tumor_id: $params.tumor_id"
println "normal_id: $params.normal_id"


Channel.fromPath( params.input_vcf, checkIfExists:true )        
	.set{ input_vcf_ch }

process vcf2maf {

  	publishDir "${params.outdir}/maf", mode:'copy'
  	
  	input:
  	file (vcf) from input_vcf_ch

  	output:
  	//set file("*.vep.vcf"), file("out_filtered_SOMATIC.maf")  into maf_ch
  	 file("${params.param}_SOMATIC.maf")  into maf_ch


  	script:
  	"""
  	echo "https://api.oncokb.org/oncokb-annotator/gene-annotation#cna-annotator"
  	echo "We recommend processing VCF files by vcf2maf with OncoKB isoforms before using the MafAnnotator here."
  	echo "We recommend processing VCF files by vcf2maf with MSK override isoforms before using the MafAnnotator here"
 	
  	#gunzip ${params.input_vcf}
   
   	perl /hpcshare/genomics/yabili_pipeline_vda/yabili_annotation_vda/vcf2maf.pl \
      	--input-vcf ${params.input_vcf} \
      	--output-maf ${params.param}_SOMATIC.maf \
      	--tumor-id ${params.tumor_id} \
      	--normal-id ${params.normal_id} \
      	--ref-fasta ${params.ref_fasta} \
      	--vep-path ${params.vep_path} \
      	--vep-data ${params.vep_data} \
      	--ncbi-build GRCh38

 	"""
  
}

process oncokb_mafannotator {

    publishDir "${params.outdir}/oncokb_mafannotator", mode:'copy'

  	input:
  	file(maf)  from maf_ch
   
  	output:
  	file("${params.param}.maf.oncokb") into oncokb_mafannotator_ch

  	script:
  	"""
  	python3 /hpcshare/genomics/yabili_pipeline_vda/yabili_annotation_vda/MafAnnotator.py \
    -i ${maf} \
    -o ${params.param}.maf.oncokb \
    -b ${TOKEN} \
    -t BRCA #in questo modo si bypassa la clinical_data
 	  """
}





/*
process clinical_data {

publishDir "${params.outdir}/...", mode:'copy'

  	input:
  	...

  	output:
  	...

  	script:
  	"""
  	Essential clinical columns:
    SAMPLE_ID: sample ID
    ONCOTREE_CODE: tumor type code from oncotree (oncotree.mskcc.org)
    Cancer type will be assigned based on the following priority:
    1) ONCOTREE_CODE in clinical data file
    2) ONCOTREE_CODE exist in MAF
    3) default tumor type (-t)

    #http://oncotree.mskcc.org/#/home?tab=mapping
 	"""
}


#https://api.oncokb.org/oncokb-annotator/gene-annotation#cna-annotator
#i tumor type che analizzeremo sono:
#BRCA
#LUNG
#COLON
#-PANCREAS
#-KIDNEY
#-MELANOMA
*/



/*
process oncokb_fusionannotator


process oncokb_cnaannotator


process oncokb_clinicaldataannotator


process oncokb_oncokbplots

*/
