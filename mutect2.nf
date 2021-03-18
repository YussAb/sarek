// Required Parameters
//https://groups.google.com/g/nextflow/c/T8l-eqc5MoI

// Print some stuff here
println "reads: $params.normal"
println "reads: $params.disease"
println "reference: $params.ref"

// Step che fa preprocessing del tsv per mutect2


//
params.tsv = "/hpcshare/genomics/yabili_pipeline_vda/yabili_germline_vda/HG002_SINGLE_GATK/results/Preprocessing/TSV/recalibrated.tsv"

Channel.fromPath( params.tsv , checkIfExists:true)
       .splitCsv(sep:'\t')
       .set { source_ch  }

normal = Channel.create()
disease = Channel.create()

source_ch.choice( normal, disease ) {row -> "${row[2]}" == "0" ? 0 : 1}

 
normal.collect { row -> "${row[4]}" }
         .into{ normal_ch ; normal_pileup_ch ; normal_pileup_bis_ch}
          
disease.collect { row -> "${row[4]}" }
         .into{ disease_ch ; disease_pileup_ch ; disease_pileup_bis_ch }             

//

 
 
 
Channel.value( params.normal_ID )        
	.into{ normal_ID_ch ; normal_ID_pileup_ch }
 
Channel.value( params.disease_ID )
	.into{ disease_ID_ch ; disease_ID_pileup_ch }


Channel.fromPath( params.pon, checkIfExists:true )
        .set{ pon_ch }

Channel.fromPath( params.pon_index, checkIfExists:true )
        .set{ pon_index_ch }


process mutect2 {

  publishDir "${params.outdir}/raw_vcf", mode: 'copy'
	
	input:
  path (normal)  from normal_ch
  //file (normal_ID)  from normal_ID_ch
	path (disease) from disease_ch
  //file (disease_ID) from disease_ID_ch
  file (pon)     from pon_ch
  file (tbi)     from pon_index_ch
	
	output:
	set file ("${params.disease_ID}_vs_${params.normal_ID}.vcf.gz"), \
	 file("${params.disease_ID}_vs_${params.normal_ID}.vcf.gz.tbi"), \
	 file("${params.disease_ID}_vs_${params.normal_ID}.vcf.gz.stats") \
	 into variant_ch
	file("f1r2.tar.gz") into f1r2_in_ch	

	script:
	"""

	gatk Mutect2 \
  -R ${params.ref} \
  -I ${disease} -tumor  ${params.disease_ID}   \
  -I ${normal}  -normal ${params.normal_ID}    \
	 --germline-resource ${params.gnomad} \
	 --panel-of-normals ${pon} --f1r2-tar-gz f1r2.tar.gz \
	 -O ${params.disease_ID}_vs_${params.normal_ID}.vcf.gz
	"""
}


process f1r2 {
	
  publishDir "${params.outdir}/f1r2", mode: 'copy'

	input:
	file(f1r2) from  f1r2_in_ch

	output:
	file ("artifact-prior.tar.gz") into f1r2_out_ch  

	script:
	"""
	gatk LearnReadOrientationModel \
	-I ${f1r2} \
	-O artifact-prior.tar.gz 
	"""
}

//https://github.com/broadinstitute/gatk-protected/blob/master/src/test/resources/org/broadinstitute/hellbender/tools/exome/common_SNP.interval_list

process getpileup_normal  {

  publishDir "${params.outdir}/contamination", mode: 'copy'

	input:
        path (normal) from normal_pileup_ch
      //file (normal_ID) from normal_ID_pileup_ch
        path (disease) from disease_pileup_ch
      //file (disease_ID) from disease_ID_pileup_ch

	output:
        file ("${params.normal_ID}_pileups.table") into normal_pileuptable_ch

	script:
	"""
	 gatk GetPileupSummaries \
	 -I ${normal} \
	 -V ${params.gnomad} \
	 -L ${params.interval} \
	 -O ${params.normal_ID}_pileups.table
  
	"""
}



process getpileup_disease {

  publishDir "${params.outdir}/contamination", mode: 'copy'

        input:
  path (normal) from normal_pileup_bis_ch
 // file (normal_ID) from normal_ID_pileup_ch
  path (disease) from disease_pileup_bis_ch
 //     file (disease_ID) from disease_ID_pileup_ch

        output:
        file ("${params.disease_ID}_pileups.table") into disease_pileuptable_ch

        script:
        """
  gatk GetPileupSummaries \
         -I ${disease} \
         -V ${params.gnomad} \
         -L ${params.interval} \
         -O ${params.disease_ID}_pileups.table
        """
}




process calculatecont {

  publishDir "${params.outdir}/contamination", mode: 'copy'

	input:
	file (disease_pileups_table) from disease_pileuptable_ch
  file (normal_pileups_table) from normal_pileuptable_ch

	output:
	file("contamination.table") into contamination_ch

	script:
	"""
 	gatk CalculateContamination \
   	-I  ${disease_pileups_table} \
    -matched ${normal_pileups_table} \
   	-O contamination.table
	"""

}


process filtervariant {

	publishDir "${params.outdir}/filt_vcf", mode: 'copy'

	input:
	set file (somatic_vcf),file(tbi),file(stats) from variant_ch
	file (contamination_table) from contamination_ch
	file (artifact_prior) from f1r2_out_ch

	output:
	file("${params.disease_ID}_vs_${params.normal_ID}_filtered.vcf.gz") into filtvcf_ch

	script:
	"""
	 gatk FilterMutectCalls \
  		-R ${params.ref} \
  		-V ${somatic_vcf} \
  		--contamination-table ${contamination_table} \
  		--ob-priors ${artifact_prior} \
  		-O ${params.disease_ID}_vs_${params.normal_ID}_filtered.vcf.gz
	"""
}

