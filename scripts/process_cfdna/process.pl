#!/usr/bin/perl -w

use strict;

use Cwd;
use FindBin;
use File::Basename;
use Data::Dumper;
use lib "$FindBin::Bin/";

use ParseArgs;
use DM;

# Take as input a MANIFEST file (a simple text file detailing some properties of the input data) and
# an output directory.
my %args = &getCommandArguments(
    'MANIFEST'    => undef,
    'OUT_DIR'     => undef,
    'DRY_RUN'     => 1,
    'NUM_JOBS'    => 1,
    'NUM_THREADS' => 1,
    'KEEP_GOING'  => 1,
    'CLEANUP'     => 1,
    'CLUSTER'     => 'localhost',
    'QUEUE'       => 'localhost',
);

# Create a job management object (this writes a Makefile to a temp file and allows the author to
# add job rules, their dependencies, and execute when finished.  It's sorta like snakemake I think.).
my $dm = new DM(
    'dryRun'     => $args{'DRY_RUN'},
    'numJobs'    => $args{'NUM_JOBS'},
    'keepGoing'  => $args{'KEEP_GOING'},
    'cluster'    => $args{'CLUSTER'},
    'queue'      => $args{'QUEUE'},
    'unlink'     => $args{'CLEANUP'},
    'outputFile' => "dm.log",
);

# Set some program paths.  We'll assume bwa, samtools, and connor-runner.py are in the PATH environment variable.
chomp(my $bwa = qx(which bwa));
chomp(my $samtools = qx(which samtools));
chomp(my $connor = qx(which connor-runner.py));
my $gatk = "$FindBin::Bin/../../gatk-launch";
my $filter = "$FindBin::Bin/filterVcf.pl";
my $ref = "$FindBin::Bin/../../resources/hs37d5.fa";
my $bed = "$FindBin::Bin/resources/my-targets-cutts.hs37d5_formatted.bed";

# We'll only take cfDNA-type projects from the MANIFEST file
my %manifest = loadManifest($args{'MANIFEST'}, 'type' => 'cfDNA');

# For each patient, we'll treat each timepoint as a separate sample, and combine the alignments in the end so
# that we can do joint calling on all of the timepoints for a patient at once.
foreach my $patient (keys(%manifest)) {
    my @entries = @{$manifest{$patient}};
    my %timepoints;

    my $rootdir = "$args{'OUT_DIR'}/$patient/";

    # First, process each sample separately.  At this stage, each sample is totally independent of the other in terms of data processing.
    foreach my $entryref (@entries) {
        my %e = %{$entryref};

        # Each hash entry from the manifest looks like this (for example):
        #{
        #  'magic' => '32726716',
        #  'project' => 'Plasma004_pre_treatment_samples-28195179',
        #  'name' => 'GEL088',
        #  'barcoded' => 'N',
        #  'fullsample' => 'GEL88-32726716',
        #  'treatment' => 'PRE',
        #  'fc' => 'S5',
        #  'type' => 'cfDNA',
        #  'end1' => '/data3/users/kiran/analyses/Schuh/data/cfDNA/Plasma004_pre_treatment_samples-28195179/GEL88-32726716/GEL88_S5_L008_R1_001.fastq.gz',
        #  'end2' => '/data3/users/kiran/analyses/Schuh/data/cfDNA/Plasma004_pre_treatment_samples-28195179/GEL88-32726716/GEL88_S5_L008_R2_001.fastq.gz',
        #  'timepoint' => '0',
        #  'lane' => '8',
        #  'material' => 'FF'
        #};

        # We'll set the read group to have enough information to definitively distinguish different runs of the same patient and timepoint.
        my $rg = "$patient.$e{'type'}.$e{'material'}.$e{'timepoint'}.$e{'fc'}.$e{'lane'}.$e{'magic'}";

        # We'll set the sample to specify the patient, the experiment type, and the timepoint.  This way, multiple runs from the same
        # experiment and timepoint will automatically register as a single sample to the GATK HaplotypeCaller, and data from different timepoints
        # will automatically register as a new sample.
        my $sm = "$patient.$e{'type'}.$e{'timepoint'}";

        # The library field in the BAM read group is unused by us, but who knows - we might use it later.
        my $lb = "$e{'type'}.$e{'material'}.$e{'fc'}.$e{'lane'}";

        # Set up the output directory and the BAM file prefix.
        my $dir = "$rootdir/$e{'timepoint'}";
        my $bamPrefix = "$dir/$rg";
        my $bam = "$bamPrefix.bam";

        # Handy one-liner to align data, set the read group, convert to BAM, and sort by coordinate.
        my $bamCmd = "$bwa mem -v 1 -t $args{'NUM_THREADS'} -R \"\@RG\\tID:$rg\\tSM:$sm\\tPL:ILLUMINA\\tLB:$lb\\tPU:$e{'fc'}.$e{'lane'}\\tDS:$e{'treatment'}\" $ref $e{'end1'} $e{'end2'} | samtools view -bS - | samtools sort - $bamPrefix";
        $dm->addRule($bam, [$e{'end1'}, $e{'end2'}], $bamCmd);

        # Index the BAM
        my $bai = "$bam.bai";
        my $baiCmd = "$samtools index $bam";
        $dm->addRule($bai, $bam, $baiCmd);

        # Run Connor on the BAM file.  Note that we're explicitly assuming from this point forward that every dataset we process will be barcoded.
        my $connorBam = "$bamPrefix.connor.bam";
        my $connorBamCmd = "$connor --force $bam $connorBam";
        $dm->addRule($connorBam, $bai, $connorBamCmd);

        # Frustratingly, Connor doesn't update the NM flag in each BAM record after emitting consensus sequences
        # to reflect the new number of mismatches.  We fix that here.
        my $calmdBam = "$bamPrefix.connor.calmd.bam";
        my $calmdBamCmd = "$samtools calmd -b $connorBam $ref > $calmdBam 2>/dev/null";
        $dm->addRule($calmdBam, $connorBam, $calmdBamCmd);

        # Index the new Connor'd BAM.
        my $calmdBai = "$calmdBam.bai";
        my $calmdBaiCmd = "$samtools index $calmdBam";
        $dm->addRule($calmdBai, $calmdBam, $calmdBaiCmd);

        push(@{$timepoints{$e{'timepoint'}}}, $calmdBam);
    }

    my @mergedTpBams;
    my @reducedTpBams;

    # Now, iterate through all the timepoints for that patient, and combine multiple datasets per timepoint into single easy-to-manipulate files.
    foreach my $timepoint (keys(%timepoints)) {
        my @tpbams = @{$timepoints{$timepoint}};

        # Merge all the data for this timepoint into a single BAM file.
        my $mergedBam = "$rootdir/$timepoint/merged.bam";
        my $mergedBamCmd = "$gatk PrintReads -R $ref -O $mergedBam -I " . join(" -I ", @tpbams);
        $dm->addRule($mergedBam, \@tpbams, $mergedBamCmd);

        # Search for reads that are the sole support for a reference mismatch in a pileup.  These are likely errors and can be removed.
        my $singleMismatchReads = "$rootdir/$timepoint/smr.txt";
        my $singleMismatchReadsCmd = "$gatk ListReadsToRemove -R $ref -O $singleMismatchReads -I " . join(" -I ", @tpbams);
        $dm->addRule($singleMismatchReads, \@tpbams, $singleMismatchReadsCmd);

        # Remove the aforementioned reads.
        my $reducedBamPrefix = "$rootdir/$timepoint/$patient.smr";
        my $reducedBam = "$reducedBamPrefix.bam";
        my $reducedBamCmd = "$gatk RemoveReads -R $ref -rl $singleMismatchReads -O $reducedBam -I " . join(" -I ", @tpbams);
        $dm->addRule($reducedBam, $singleMismatchReads, $reducedBamCmd);

        # Now remove reads that are a perfect match to the reference.
        my $moreReducedBamPrefix = "$reducedBamPrefix.pr";
        my $moreReducedBam = "$moreReducedBamPrefix.bam";
        my $moreReducedBamCmd = "$samtools view -h $reducedBam | grep -v 'NM:i:0' | samtools view -bS - | samtools sort - $moreReducedBamPrefix && samtools index $moreReducedBam";
        $dm->addRule($moreReducedBam, $reducedBam, $moreReducedBamCmd);

        # These are the tumor+normal cfDNA BAM files
        push(@mergedTpBams, $mergedBam);

        # These are the tumor-only cfDNA BAM files (provided we've done the germline subtraction properly)
        push(@reducedTpBams, $moreReducedBam);
    }

    # Merge the tumor-only cfDNA data from all timepoints into a single file.
    my $reducedFinal = "$rootdir/$patient.glsubtract.bam";
    my $reducedFinalCmd = "$gatk PrintReads -R $ref -O $reducedFinal -I " . join(" -I ", sort { $a cmp $b } @reducedTpBams);
    $dm->addRule($reducedFinal, \@reducedTpBams, $reducedFinalCmd);

    # Merge the tumor+normal cfDNA data from all timepoints into a single file.
    my $mergedFinal = "$rootdir/$patient.full.bam";
    my $mergedFinalCmd = "$gatk PrintReads -R $ref -O $mergedFinal -I " . join(" -I ", sort { $a cmp $b } @mergedTpBams);
    $dm->addRule($mergedFinal, \@mergedTpBams, $mergedFinalCmd);

    # Call the GATK's haplotype caller on the tumor-only data.
    my $reducedVcf = "$rootdir/$patient.glsubtract.vcf";
    my $reducedVcfCmd = "$gatk HaplotypeCaller -R $ref -I $reducedFinal -L $bed -O $reducedVcf --dontUseSoftClippedBases --maxReadsPerAlignmentStart 1000 -errorCorrectReads --min_base_quality_score 25 -RF AmbiguousBaseReadFilter -ploidy 1 -threads 12";
    $dm->addRule($reducedVcf, $reducedFinal, $reducedVcfCmd);

    # Annotate the tumor-only VCF with the read depth data from the ***TUMOR+NORMAL*** data.
    # This step serves two purposes:
    # 1. later, when we try to compute the allele fraction of the putative tumor allele, we need the depth
    #    in terms of the unsubtracted data, but because we ran HaplotypeCaller on the subtracted data, this
    #    needs to be added back as a post-processing step.
    # 2. there is a bug in GATK4 wherein VCF indices are not currently emitted by HaplotypeCaller, which breaks
    #    usage of the VCF with other GATK tools.  Running the tool below has the side effect of creating an index.
    my $reducedAnnotatedVcf = "$rootdir/$patient.glsubtract.annotated.vcf";
    my $reducedAnnotatedVcfCmd = "$gatk AnnotateVcfWithBamDepth -R $ref -I $mergedFinal -V $reducedVcf -O $reducedAnnotatedVcf";
    $dm->addRule($reducedAnnotatedVcf, $reducedVcf, $reducedAnnotatedVcfCmd);

    # Call the GATK's haplotype caller on the tumor+normal data.  Here, we make the assumption that the ctDNA fraction is so low, it
    # will look like sequencing error to the caller and thus be ignored.  This will give us a VCF of het or hom-var sites in the patient
    # that we can ignore.
    my $mergedVcf = "$rootdir/$patient.full.vcf";
    my $mergedVcfCmd = "$gatk HaplotypeCaller -R $ref -I $mergedFinal -L $reducedVcf --maxReadsPerAlignmentStart 1000 -activeProbabilityThreshold 0 -stand_call_conf 1 -threads 12 -O $mergedVcf";
    $dm->addRule($mergedVcf, $mergedFinal, $mergedVcfCmd);

    # Annotate BAM depth.
    my $mergedAnnotatedVcf = "$rootdir/$patient.full.annotated.vcf";
    my $mergedAnnotatedVcfCmd = "$gatk AnnotateVcfWithBamDepth -R $ref -I $mergedFinal -V $mergedVcf -O $mergedAnnotatedVcf";
    $dm->addRule($mergedAnnotatedVcf, $mergedVcf, $mergedAnnotatedVcfCmd);

    # Compute the disjoint between the tumor-only callset and the tumor-normal callset, and retain only those sites
    # unique to the tumor.
    my $uniqueVcf = "$rootdir/$patient.unique.vcf";
    my $uniqueVcfCmd = "$gatk SelectVariants -R $ref -V $reducedAnnotatedVcf -disc $mergedAnnotatedVcf -O $uniqueVcf";
    $dm->addRule($uniqueVcf, [$reducedAnnotatedVcf, $mergedAnnotatedVcf], $uniqueVcfCmd);

    # Filter out many calls.  This is based on a very simple heuristic: does the allele fraction grow as a function of timepoint.
    # Surprisingly, this seems to remove a lot of cruft from the calls, but we need to examine this step much more closely.  There
    # are far better ways to achieve what I'm trying to do here.
    #
    # Note: this step will almost certainly break for single-timepoint data.
    my $filterVcf = "$rootdir/$patient.unique.filtered.vcf";
    my $filterVcfCmd = "$filter VCF=$uniqueVcf OUT=$filterVcf";
    $dm->addRule($filterVcf, $uniqueVcf, $filterVcfCmd);

    # We should annotate the VCF with variant presence in COSMIC and other data sources here.  The GATK tool to do this is currently broken.
    # We should switch to vcftools or something like that instead, as the GATK no longer seems to be the go-to toolkit for manipulating
    # VCF files in this manner.
    #my $cosmic = "$rootdir/$patient.unique.filtered.cosmic.vcf";
    #my $cosmicCmd = "$gatk CombineVariants -R $ref -V $filterVcf -V resources/CosmicCodingMuts.vcf -L $bed -o $cosmic";
    #$dm->addRule($cosmic, $filterVcf, $cosmicCmd);
}

$dm->execute();

# This function loads a manifest file and optionally subselects rows within it.
sub loadManifest {
    my ($manifest, %filters) = @_;

    open(MANIFEST, $manifest);
    chomp(my $header = <MANIFEST>);
    my @header = split(/\s+/, $header);

    my %manifest;

    while (my $line = <MANIFEST>) {
        chomp($line);
        my @field = split(/\s+/, $line);

        my %e;
        my $keep = 1;
        for (my $i = 0; $i <= $#header; $i++) {
            $e{$header[$i]} = $field[$i];

            if (scalar(@_) == 1 || (exists($filters{$header[$i]}) && $field[$i] !~ /$filters{$header[$i]}/)) {
                $keep = 0;
            }
        }

        if ($keep) {
            push(@{$manifest{$e{'name'}}}, \%e);
        }
    }
    close(MANIFEST);

    return %manifest;
}
