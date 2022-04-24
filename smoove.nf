#!/usr/bin/env nextflow

   /* CALL CNVS USING SMOOVE (modified lumpy for groups)
   ================================================================================
           the required command line inputs are:
                 run      nf file
                 -c       config file
                 --ref    reference genome file
                 --list   file of BAMs formatted: SAMPLEID, FILE DIRECTORY
//               --chr    file lising chromosome names
                 --out     directory for publishing output
  ===============================================================================
 */
t
//      INITIALIZE INPUT

DAT = file("${params.list}") // #input should <SAMPLEID>,<BAM FILE DIRECTORY>,
REF = file("${params.ref}")
OUT_DIR = "${params.out}"

smoove_SING = "singularity exec /beegfs/work/public/singularity/smoove-0.2.3.simg"

Channel
        .fromPath("${DAT}")
        .splitCsv(header: ['SAMP_ID', 'FILE'])
        .map{ row -> tuple(row.SAMP_ID, row.FILE, file(row.FILE)) }
        .set{sampledata}

// =========== call genotypes for each sample with smoove ==============

process callgeno {
"CNV_nextflow_smoove.nf" 126L, 3738C                                             1,1           Top
#!/usr/bin/env nextflow

   /* CALL CNVS USING SMOOVE (modified lumpy for groups)
   ================================================================================

           the required command line inputs are:
                 run      nf file
                 -c       config file
                 --ref    reference genome file
                 --list   file of BAMs formatted: SAMPLEID, FILE DIRECTORY
//               --chr    file lising chromosome names
                 --out     directory for publishing output
  ===============================================================================
 */
t
//      INITIALIZE INPUT
REF = file("${params.ref}")
OUT_DIR = "${params.out}"

Channel
        .fromPath("${DAT}")
        .splitCsv(header: ['SAMP_ID', 'FILE'])
        .map{ row -> tuple(row.SAMP_ID, row.FILE, file(row.FILE)) }
        .set{sampledata}

// =========== call genotypes for each sample with smoove ==============

process callgeno {

        tag {ID}
        errorStrategy 'retry'
        maxRetries 2
        executor = 'slurm'

        clusterOptions = "--cpus-per-task=${params.cpu} --time=${params.time} --mem=${params.mem}"

        input:
        set val(ID), val(BAM_path), file(BAM_FILE) from sampledata

        output:
        file({ "${ID}-smoove.genotyped.vcf.gz" }) into GENOCALL
        set file({ "${BAM_FILE}" }), val(BAM_path), val(ID) into CALL_GENO_OUT

        script:
        """
        ${smoove_SING} smoove call --outdir . --name $ID --fasta $REF -p 1 --genotype ${BAM_path}
        """
}

         errorStrategy 'retry'
         maxRetries 2

         input:
         file allgeno from GENOCALL.collect()

         output:
         file( {"unified-samples.sites.vcf.gz"} ) into SMOOVE_UNIFIED
         """
         ${smoove_SING} smoove merge --name unified-samples --fasta ${REF} --outdir ./ ./${allgeno}
         """

UNIFIED= CALL_GENO_OUT.combine(SMOOVE_UNIFIED)

// ======== genotype each sample from smoove (makes VCF files) ======
// smoove inputs: name, *test-mergedsites.vcf BAM file, makes VCF for each sample

process genotype_smoove {

        tag { ID }
        errorStrategy 'retry'
        maxRetries 2
        executor = 'slurm'

        input:
        set file(BAM_FILE), val(BAM_PATH), val(ID), file(UNIFIED_VCF) from UNIFIED

        output:
        set val(ID), file({ "${BAM_FILE}" }) into CNVnator_INPUT
        file("$ID-joint-smoove.genotyped.vcf.gz") into SMOOVE_OUT
        file("$ID-joint-smoove.genotyped.vcf.gz.csi") into CZI_files

        script:
        """
        ${smoove_SING} smoove genotype -x -d -p 1 --name ${ID}-joint --outdir ./ --fasta ${REF} --vcf ${UNIFIED_VCF} ${BAM_PATH}
        """
}

process paste_smoove {

        tag { paste_smoove }
        errorStrategy 'retry'
        maxRetries 2
        executor = 'slurm'
        clusterOptions = "--cpus-per-task=${params.cpu} --time=${params.time} --mem=${params.mem}"

        publishDir "${OUT_DIR}", mode: 'move'

        input:
        file (allgeno) from SMOOVE_OUT.collect()
        file(allczi) from CZI_files.collect()

        output:
        file ("smoove_out*") into SMOOVE_VCF

        script:
        println "geno is ${allgeno}"
        println "czi is ${allczi}"
        """
        ${smoove_SING} smoove paste --name smoove_out ${allgeno}
        """
}
