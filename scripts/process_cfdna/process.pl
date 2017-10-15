#!/usr/bin/perl -w

use strict;

use Cwd;
use FindBin;
use File::Basename;
use Data::Dumper;
use lib "$FindBin::Bin/";

use ParseArgs;
use DM;

my %args = &getCommandArguments(
    'MANIFEST'    => undef,
    'OUT_DIR'     => undef,
    'REF'         => undef,
    'DRY_RUN'     => 1,
    'NUM_JOBS'    => 1,
    'NUM_THREADS' => 1,
    'KEEP_GOING'  => 1,
    'CLUSTER'     => 'localhost',
    'QUEUE'       => 'localhost',
);

my $dm = new DM(
    'dryRun'     => $args{'DRY_RUN'},
    'numJobs'    => $args{'NUM_JOBS'},
    'keepGoing'  => $args{'KEEP_GOING'},
    'cluster'    => $args{'CLUSTER'},
    'queue'      => $args{'QUEUE'},
    'outputFile' => "dm.log",
);

my $bwa = "/home/kiran/repositories/bwa/bwa";
my $samtools = "/home/kiran/repositories/samtools/samtools";
my $connor = "/home/kiran/repositories/Connor/connor-runner.py";
my $gatk = "/home/kiran/repositories/CtDNATools/gatk-launch";
my $filter = "/home/kiran/repositories/CtDNATools/scripts/process_cfdna/filterVcf.pl";

my %manifest = loadManifest($args{'MANIFEST'}, 'type' => 'cfDNA', 'name' => 'GEL212');

foreach my $patient (keys(%manifest)) {
    my @entries = @{$manifest{$patient}};
    my %timepoints;

    my $rootdir = "$args{'OUT_DIR'}/$patient/";

    foreach my $entryref (@entries) {
        my %e = %{$entryref};

        my $rg = "$patient.$e{'type'}.$e{'material'}.$e{'timepoint'}.$e{'fc'}.$e{'lane'}.$e{'magic'}";
        my $sm = "$patient.$e{'type'}.$e{'timepoint'}";
        my $lb = "$e{'type'}.$e{'material'}.$e{'fc'}.$e{'lane'}";

        my $oldBam = "results/process/$patient/$e{'type'}/$e{'treatment'}/$e{'timepoint'}/$patient.$e{'fc'}.$e{'lane'}.$e{'barcoded'}.$e{'magic'}.bam";

        my $dir = "$rootdir/$e{'timepoint'}";
        my $bamPrefix = "$dir/$rg";
        my $bam = "$bamPrefix.bam";

        if (-e $oldBam) {
            my $header = "$bamPrefix.header";
            my $headerCmd = "samtools view -H $oldBam | sed 's/^\@RG.\\+/\@RG\\tID:$rg\\tSM:$sm\\tPL:ILLUMINA\\tLB:$lb\\tPU:$e{'fc'}.$e{'lane'}\\tDS:$e{'treatment'}/' > $header";
            $dm->addRule($header, $oldBam, $headerCmd);

            my $bamCmd = "samtools reheader $header $oldBam > $bam";
            $dm->addRule($bam, $header, $bamCmd);
        } else {
            my $bamCmd = "$bwa mem -v 1 -t $args{'NUM_THREADS'} -R \"\@RG\\tID:$rg\\tSM:$sm\\tPL:ILLUMINA\\tLB:$lb\\tPU:$e{'fc'}.$e{'lane'}\\tDS:$e{'treatment'}\" $args{'REF'} $e{'end1'} $e{'end2'} | samtools view -bS - | samtools sort - $bamPrefix";
            $dm->addRule($bam, [$e{'end1'}, $e{'end2'}], $bamCmd);
        }

        my $bai = "$bam.bai";
        my $baiCmd = "$samtools index $bam";
        $dm->addRule($bai, $bam, $baiCmd);

        my $connorBam = "$bamPrefix.connor.bam";
        my $connorBamCmd = "$connor --force $bam $connorBam";
        $dm->addRule($connorBam, $bai, $connorBamCmd);

        my $calmdBam = "$bamPrefix.connor.calmd.bam";
        my $calmdBamCmd = "$samtools calmd -b $connorBam $args{'REF'} > $calmdBam 2>/dev/null";
        $dm->addRule($calmdBam, $connorBam, $calmdBamCmd);

        my $calmdBai = "$calmdBam.bai";
        my $calmdBaiCmd = "$samtools index $calmdBam";
        $dm->addRule($calmdBai, $calmdBam, $calmdBaiCmd);

        my $repaired = "$bamPrefix.connor.calmd.repaired.bam";
        my $repairedCmd = "$samtools addreplacerg -m overwrite_all -R $rg -o $repaired $calmdBam";
        $dm->addRule($repaired, $calmdBai, $repairedCmd);

        push(@{$timepoints{$e{'timepoint'}}}, $repaired);
    }

    my @mergedTpBams;
    my @reducedTpBams;

    foreach my $timepoint (keys(%timepoints)) {
        my @tpbams = @{$timepoints{$timepoint}};

        my $mergedBam = "$rootdir/$timepoint/merged.bam";
        my $mergedBamCmd = "$gatk PrintReads -R $args{'REF'} -O $mergedBam -I " . join(" -I ", @tpbams);
        $dm->addRule($mergedBam, \@tpbams, $mergedBamCmd);

        my $singleMismatchReads = "$rootdir/$timepoint/smr.txt";
        my $singleMismatchReadsCmd = "$gatk ListReadsToRemove -R $args{'REF'} -O $singleMismatchReads -I " . join(" -I ", @tpbams);
        $dm->addRule($singleMismatchReads, \@tpbams, $singleMismatchReadsCmd);

        my $reducedBamPrefix = "$rootdir/$timepoint/$patient.smr";
        my $reducedBam = "$reducedBamPrefix.bam";
        my $reducedBamCmd = "$gatk RemoveReads -R $args{'REF'} -rl $singleMismatchReads -O $reducedBam -I " . join(" -I ", @tpbams);
        $dm->addRule($reducedBam, $singleMismatchReads, $reducedBamCmd);

        my $moreReducedBamPrefix = "$reducedBamPrefix.pr";
        my $moreReducedBam = "$moreReducedBamPrefix.bam";
        my $moreReducedBamCmd = "$samtools view -h $reducedBam | grep -v 'NM:i:0' | samtools view -bS - | samtools sort - $moreReducedBamPrefix && samtools index $moreReducedBam";
        $dm->addRule($moreReducedBam, $reducedBam, $moreReducedBamCmd);

        push(@mergedTpBams, $mergedBam);
        push(@reducedTpBams, $moreReducedBam);
    }

    my $reducedFinal = "$rootdir/$patient.glsubtract.bam";
    my $reducedFinalCmd = "$gatk PrintReads -R $args{'REF'} -O $reducedFinal -I " . join(" -I ", sort { $a cmp $b } @reducedTpBams);
    $dm->addRule($reducedFinal, \@reducedTpBams, $reducedFinalCmd);

    my $mergedFinal = "$rootdir/$patient.full.bam";
    my $mergedFinalCmd = "$gatk PrintReads -R $args{'REF'} -O $mergedFinal -I " . join(" -I ", sort { $a cmp $b } @mergedTpBams);
    $dm->addRule($mergedFinal, \@mergedTpBams, $mergedFinalCmd);

    my $reducedVcf = "$rootdir/$patient.glsubtract.vcf";
    my $reducedVcfCmd = "$gatk HaplotypeCaller -R $args{'REF'} -I $reducedFinal -L resources/my-targets-cutts.hs37d5_formatted.bed -O $reducedVcf --dontUseSoftClippedBases --maxReadsPerAlignmentStart 1000 -errorCorrectReads --min_base_quality_score 25 -RF AmbiguousBaseReadFilter -ploidy 1 -threads 12";
    $dm->addRule($reducedVcf, $reducedFinal, $reducedVcfCmd);

    my $reducedAnnotatedVcf = "$rootdir/$patient.glsubtract.annotated.vcf";
    my $reducedAnnotatedVcfCmd = "$gatk AnnotateVcfWithBamDepth -R $args{'REF'} -I $mergedFinal -V $reducedVcf -O $reducedAnnotatedVcf";
    $dm->addRule($reducedAnnotatedVcf, $reducedVcf, $reducedAnnotatedVcfCmd);

    my $mergedVcf = "$rootdir/$patient.full.vcf";
    my $mergedVcfCmd = "$gatk HaplotypeCaller -R $args{'REF'} -I $mergedFinal -L $reducedVcf --maxReadsPerAlignmentStart 1000 -activeProbabilityThreshold 0 -stand_call_conf 1 -threads 12 -O $mergedVcf";
    $dm->addRule($mergedVcf, $mergedFinal, $mergedVcfCmd);

    my $mergedAnnotatedVcf = "$rootdir/$patient.full.annotated.vcf";
    my $mergedAnnotatedVcfCmd = "$gatk AnnotateVcfWithBamDepth -R $args{'REF'} -I $mergedFinal -V $mergedVcf -O $mergedAnnotatedVcf";
    $dm->addRule($mergedAnnotatedVcf, $mergedVcf, $mergedAnnotatedVcfCmd);

    my $uniqueVcf = "$rootdir/$patient.unique.vcf";
    my $uniqueVcfCmd = "$gatk SelectVariants -R $args{'REF'} -V $reducedAnnotatedVcf -disc $mergedAnnotatedVcf -O $uniqueVcf";
    $dm->addRule($uniqueVcf, [$reducedAnnotatedVcf, $mergedAnnotatedVcf], $uniqueVcfCmd);

    my $filterVcf = "$rootdir/$patient.unique.filtered.vcf";
    my $filterVcfCmd = "$filter VCF=$uniqueVcf OUT=$filterVcf";
    $dm->addRule($filterVcf, $uniqueVcf, $filterVcfCmd);

    my $cosmic = "$rootdir/$patient.unique.filtered.cosmic.vcf";
    my $cosmicCmd = "$gatk CombineVariants -R $args{'REF'} -V $filterVcf -V resources/CosmicCodingMuts.vcf -L resources/my-targets-cutts.hs37d5_formatted.bed -o $cosmic";
    $dm->addRule($cosmic, $filterVcf, $cosmicCmd);
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
