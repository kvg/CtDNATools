package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import avro.shaded.com.google.common.annotations.VisibleForTesting;
import com.google.cloud.genomics.dataflow.utils.GCSOptions;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.RDDUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SerializablePredicate;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;
import scala.Tuple3;

import java.util.*;
import java.util.stream.Collectors;


final class ForSimpleStrandSwitch implements VariantDetectorFromLongReadAlignments {

    @SuppressWarnings("unchecked")
    private static final List<String> EMPTY_INSERTION_MAPPINGS = Collections.EMPTY_LIST;

    static final int MORE_RELAXED_ALIGNMENT_MIN_LENGTH = 30;
    static final int MORE_RELAXED_ALIGNMENT_MIN_MQ = 20;

    public void inferSvAndWriteVCF(final JavaRDD<AlignedContig> longReads, final String vcfOutputFileName,
                                   final Broadcast<ReferenceMultiSource> broadcastReference, final String fastaReference,
                                   final GCSOptions options, final Logger toolLogger) {

        longReads.cache();
        toolLogger.info(longReads.count() + " chimera indicating either 1) simple strand-switch breakpoints, or 2) inverted duplication.");

        // split between suspected inv dup VS strand-switch breakpoint
        final Tuple2<JavaRDD<AlignedContig>, JavaRDD<AlignedContig>> split
                = RDDUtils.split(longReads.map(ForSimpleStrandSwitch::removeOverlap),
                                  new IsLikelyInvertedDuplication(), true);

        final JavaRDD<VariantContext> simpleStrandSwitchBkpts =
                dealWithSimpleStrandSwitchBkpts(split._2, broadcastReference, toolLogger);
        SVVCFWriter.writeVCF(options, vcfOutputFileName.replace(".vcf", "_simpleSS.vcf"),
                fastaReference, simpleStrandSwitchBkpts, toolLogger);

//        final JavaRDD<VariantContext> invDups =
//                dealWithSuspectedInvDup(split._1, broadcastReference, toolLogger);
//        if (invDups != null)
//            SVVCFWriter.writeVCF(options, vcfOutputFileName.replace(".vcf", "_invDup.vcf"),
//                    fastaReference, invDups, toolLogger);
    }

    public static final class IsLikelyInvertedDuplication implements SerializablePredicate<AlignedContig> {
        private static final long serialVersionUID = 1L;
        @Override
        public boolean test(final AlignedContig longRead) {
            return BreakpointComplications.isLikelyInvertedDuplication(longRead.alignmentIntervals.get(0),
                                                                       longRead.alignmentIntervals.get(1));
        }
    }

    // =================================================================================================================

    /**
     * Removes overlap from a designated alignment interval, so that the inverted duplicated reference span is minimal.
     * If the two alignment intervals are NOT overlapping, return the original read.
     */
    private static AlignedContig removeOverlap(final AlignedContig longRead) {
        final int overlapOnRead = AlignmentInterval.overlapOnContig(longRead.alignmentIntervals.get(0),
                                                                    longRead.alignmentIntervals.get(1));
        if (overlapOnRead==0) {
            return longRead;
        } else {
            final AlignmentInterval one = longRead.alignmentIntervals.get(0),
                                    two = longRead.alignmentIntervals.get(1);
            final int js = one.referenceSpan.getEnd(),
                      jl = two.referenceSpan.getStart();
            final AlignmentInterval reconstructedOne, reconstructedTwo;
            if (js <= jl ^ one.forwardStrand) {
                reconstructedOne = one;
                reconstructedTwo = clipAlignmentInterval(two, overlapOnRead, false);
            } else {
                reconstructedOne = clipAlignmentInterval(one, overlapOnRead, true);
                reconstructedTwo = two;
            }
            return new AlignedContig(longRead.contigName, longRead.contigSequence,
                                     Arrays.asList(reconstructedOne, reconstructedTwo),
                                     longRead.hasEquallyGoodAlnConfigurations);
        }
    }

    private static AlignmentInterval clipAlignmentInterval(final AlignmentInterval input, final int clipLength,
                                                           final boolean clipFrom3PrimeEnd) {
        Utils.validateArg(clipLength < input.endInAssembledContig - input.startInAssembledContig + 1,
                            "input alignment to be clipped away: " + input.toPackedString() + "\twith clip length: " + clipLength);

        final Tuple2<SimpleInterval, Cigar> result = computeNewRefSpanAndCigar(input, clipLength, clipFrom3PrimeEnd);
        return new AlignmentInterval(result._1, input.startInAssembledContig, input.endInAssembledContig - clipLength,
                result._2, input.forwardStrand, input.mapQual,
                StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.MISSING_NM,
                input.alnScore, input.isFromSplitGapAlignment, input.hasUndergoneOverlapRemoval);
    }

    @VisibleForTesting
    static Tuple2<SimpleInterval, Cigar> computeNewRefSpanAndCigar(final AlignmentInterval input, final int clipLength,
                                                                   final boolean clipFrom3PrimeEnd) {
        Utils.validateArg(input.cigarAlong5to3DirectionOfContig.getCigarElements().stream().map(CigarElement::getOperator)
                        .noneMatch(op -> op.equals(CigarOperator.N) || op.isPadding()),
                "Input alignment contains padding or skip operations, which is currently unsupported: " + input.toPackedString());

        final Tuple3<List<CigarElement>, List<CigarElement>, List<CigarElement>> threeSections = extractCigarElements(input);

        final List<CigarElement> cigarElements = threeSections._2();
        int idx;
        final int step, stop;
        if (clipFrom3PrimeEnd) {
            idx = cigarElements.size() - 1; step = -1; stop = -1;
        } else {
            idx = 0; step = 1; stop = cigarElements.size();
        }

        Cigar newCigar = new Cigar();
        int readBasesConsumed = 0, refBasesConsumed = 0;
        while (stop != idx) {
            final CigarElement ce = cigarElements.get(idx);
            if ( ce.getOperator().equals(CigarOperator.HARD_CLIP))
                continue;
            if ( ce.getOperator().consumesReadBases() ) {
                if (readBasesConsumed + ce.getLength() < clipLength)
                    readBasesConsumed += ce.getLength();
                else { // enough read bases would be clipped

                    if ( !ce.getOperator().isAlignment() && !ce.getOperator().equals(CigarOperator.I))
                        throw new GATKException.ShouldNeverReachHereException("Logic error, should not reach here");

                    // dead with cigar first
                    final List<CigarElement> resultCEs = threeSections._1();
                    final int a = readBasesConsumed + ce.getLength() - clipLength;
                    final CigarOperator op = ce.getOperator().isAlignment() ? CigarOperator.M : CigarOperator.S;
                    if (clipFrom3PrimeEnd) {
                        resultCEs.addAll(cigarElements.subList(0, idx));
                        if (a!=0) {
                            resultCEs.add( new CigarElement(a, op) );
                        }
                        resultCEs.add(new CigarElement(clipLength, CigarOperator.S));
                    } else {
                        resultCEs.add(new CigarElement(clipLength, CigarOperator.S));
                        if (a!=0) {
                            resultCEs.add( new CigarElement(a, op) );
                        }
                        resultCEs.addAll(cigarElements.subList(idx+1, cigarElements.size()));
                    }
                    if (!threeSections._3().isEmpty())
                        resultCEs.addAll(threeSections._3());
                    newCigar = new Cigar(SvCigarUtils.compactifyNeighboringSoftClippings(resultCEs));

                    // then deal with ref span, note that here we can have only either 'M' or 'I'
                    refBasesConsumed += ce.getOperator().isAlignment() ? (clipLength - readBasesConsumed)
                                                                       : ce.getLength();

                    break;
                }
            }
            if ( ce.getOperator().consumesReferenceBases() ) { // if reaches here, not enough read bases have been consumed
                refBasesConsumed += ce.getLength();
            }
            idx += step;
        }

        if (newCigar.getCigarElements().isEmpty())
            throw new GATKException("Logic error: new cigar is empty.\t" + input.toPackedString() + "\tclip length " +
                    clipLength + "\tclip from end " + (clipFrom3PrimeEnd? "3":"5"));

        final SimpleInterval newRefSpan;
        if (clipFrom3PrimeEnd == input.forwardStrand) {
            newRefSpan = new SimpleInterval(input.referenceSpan.getContig(), input.referenceSpan.getStart(),
                    input.referenceSpan.getEnd() - refBasesConsumed);
        } else {
            newRefSpan = new SimpleInterval(input.referenceSpan.getContig(), input.referenceSpan.getStart() + refBasesConsumed,
                    input.referenceSpan.getEnd());
        }

        return new Tuple2<>(newRefSpan, newCigar);
    }

    private static Tuple3<List<CigarElement>, List<CigarElement>, List<CigarElement>> extractCigarElements(final AlignmentInterval input) {
        final List<CigarElement> cigarElements = input.cigarAlong5to3DirectionOfContig.getCigarElements();
        final List<CigarElement> left = new ArrayList<>(cigarElements.size()),
                                 middle = new ArrayList<>(cigarElements.size()),
                                 right = new ArrayList<>(cigarElements.size());
        int readBasesConsumed = 0;
        for(final CigarElement ce : cigarElements) {
            if (readBasesConsumed < input.startInAssembledContig-1) {
                left.add(ce);
            } else if (readBasesConsumed < input.endInAssembledContig) {
                middle.add(ce);
            } else {
                right.add(ce);
            }
            readBasesConsumed += ce.getOperator().consumesReadBases() || ce.getOperator().equals(CigarOperator.H)? ce.getLength() : 0;
        }

        if (middle.isEmpty())
            throw new GATKException("Logic error: cigar elements corresponding to alignment block is empty. " + input.toPackedString());

        return new Tuple3<>(left, middle, right);
    }

    /**
     * Taking advantage of the fact that for input read, we know it has only two alignments that map to the same reference
     * chromosome, with strand switch.
     * @return  null if the pair of the alignments are no strong enough to support a strand switch breakpoint
     *                  {@link #splitPairStrongEnoughEvidenceForCA(AlignmentInterval, AlignmentInterval, int, int)},
     *          otherwise a pair {chimeric alignment, read sequence}
     */
    private static Tuple2<ChimericAlignment, byte[]> convertAlignmentIntervalToChimericAlignment
    (final AlignedContig longReadWith2AIMappedToSameChrAndStrandSwitch) {

        final AlignmentInterval intervalOne = longReadWith2AIMappedToSameChrAndStrandSwitch.alignmentIntervals.get(0),
                intervalTwo = longReadWith2AIMappedToSameChrAndStrandSwitch.alignmentIntervals.get(1);

        if (splitPairStrongEnoughEvidenceForCA(intervalOne, intervalTwo, MORE_RELAXED_ALIGNMENT_MIN_MQ,  MORE_RELAXED_ALIGNMENT_MIN_LENGTH)) {
            return new Tuple2<>(new ChimericAlignment(intervalOne, intervalTwo, EMPTY_INSERTION_MAPPINGS,
                    longReadWith2AIMappedToSameChrAndStrandSwitch.contigName), longReadWith2AIMappedToSameChrAndStrandSwitch.contigSequence);
        } else {
            return null;
        }
    }

    /**
     * Roughly similar to {@link ChimericAlignment#nextAlignmentMayBeNovelInsertion(AlignmentInterval, AlignmentInterval, Integer)}:
     *  1) either alignment may have very low mapping quality (a more relaxed mapping quality threshold);
     *  2) either alignment may consume only a "short" part of the contig, or if assuming that the alignment consumes
     *     roughly the same amount of ref bases and read bases, has isAlignment that is too short
     */
    private static boolean splitPairStrongEnoughEvidenceForCA(final AlignmentInterval intervalOne,
                                                              final AlignmentInterval intervalTwo,
                                                              final int mapQThresholdInclusive,
                                                              final int alignmentLengthThresholdInclusive) {

        if (intervalOne.mapQual < mapQThresholdInclusive || intervalTwo.mapQual < mapQThresholdInclusive)
            return false;

        final int overlap = AlignmentInterval.overlapOnContig(intervalOne, intervalTwo);

        final int x = intervalOne.endInAssembledContig - intervalOne.startInAssembledContig + 1,
                y = intervalTwo.endInAssembledContig - intervalTwo.startInAssembledContig + 1;

        return Math.min(x - overlap, y - overlap) >= alignmentLengthThresholdInclusive;
    }

    // =================================================================================================================

    private JavaRDD<VariantContext> dealWithSimpleStrandSwitchBkpts(final JavaRDD<AlignedContig> longReads,
                                                                    final Broadcast<ReferenceMultiSource> broadcastReference,
                                                                    final Logger toolLogger) {

        final JavaPairRDD<ChimericAlignment, byte[]> simpleStrandSwitchBkpts =
                longReads
                        .mapToPair(ForSimpleStrandSwitch::convertAlignmentIntervalToChimericAlignment)
                        .filter(Objects::nonNull).cache();

        toolLogger.info(simpleStrandSwitchBkpts.count() + " chimera indicating simple strand-switch breakpoints.");

        return simpleStrandSwitchBkpts
                .mapToPair(pair -> new Tuple2<>(new NovelAdjacencyReferenceLocations(pair._1, pair._2), pair._1))
                .groupByKey()
                .mapToPair(noveltyAndEvidence -> inferBNDType(noveltyAndEvidence, broadcastReference.getValue()))
                .flatMap(noveltyTypeAndEvidence ->
                        AnnotatedVariantProducer
                                .produceMultipleAnnotatedVcFromNovelAdjacency(noveltyTypeAndEvidence._1,
                                        noveltyTypeAndEvidence._2._1, noveltyTypeAndEvidence._2._2, broadcastReference));
    }

    private static Tuple2<NovelAdjacencyReferenceLocations, Tuple2<Iterable<SvType>, Iterable<ChimericAlignment>>>
    inferBNDType(final Tuple2<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> noveltyAndEvidence,
                 final ReferenceMultiSource reference) {

        final NovelAdjacencyReferenceLocations novelAdjacency = noveltyAndEvidence._1;
        final Iterable<ChimericAlignment> chimericAlignments = noveltyAndEvidence._2;
        final BreakEndVariantType bkpt_1, bkpt_2;
        if (novelAdjacency.strandSwitch == StrandSwitch.FORWARD_TO_REVERSE) {
            bkpt_1 = new BreakEndVariantType.INV55BND(novelAdjacency, true, reference);
            bkpt_2 = new BreakEndVariantType.INV55BND(novelAdjacency, false, reference);
        } else if (novelAdjacency.strandSwitch == StrandSwitch.REVERSE_TO_FORWARD){
            bkpt_1 = new BreakEndVariantType.INV33BND(novelAdjacency, true, reference);
            bkpt_2 = new BreakEndVariantType.INV33BND(novelAdjacency, false, reference);
        } else {
            throw new GATKException("Wrong type of novel adjacency sent to wrong analysis pathway: no strand-switch being sent to strand-switch path. \n" +
                    Utils.stream(chimericAlignments).map(ChimericAlignment::onErrStringRep).collect(Collectors.toList()));
        }

        return new Tuple2<>(novelAdjacency, new Tuple2<>(Arrays.asList(bkpt_1, bkpt_2), chimericAlignments));
    }

    // =================================================================================================================

//    private JavaRDD<VariantContext> dealWithSuspectedInvDup(final JavaRDD<AlignedContig> longReads,
//                                                            final Broadcast<ReferenceMultiSource> broadcastReference,
//                                                            final Logger toolLogger) {
//
//        final JavaPairRDD<ChimericAlignment, byte[]> invDupSuspects =
//                longReads
//                        .mapToPair(ForSimpleStrandSwitch::convertAlignmentIntervalToChimericAlignment)
//                        .filter(Objects::nonNull).cache();
//
//        toolLogger.info(invDupSuspects.count() + " chimera indicating inverted duplication");
//
//        return invDupSuspects
//                .mapToPair(pair -> new Tuple2<>(new NovelAdjacencyReferenceLocations(pair._1, pair._2), pair._1))
//                .groupByKey()
//                .mapToPair(noveltyAndEvidence -> inferInvDupRange(noveltyAndEvidence, broadcastReference.getValue()))
//                .map(noveltyTypeAndEvidence ->
//                        DiscoverVariantsFromContigAlignmentsSAMSpark
//                                .annotateVariant(noveltyTypeAndEvidence._1, noveltyTypeAndEvidence._2._1,
//                                                 noveltyTypeAndEvidence._2._2, broadcastReference));
//    }
//
//    private static Tuple2<NovelAdjacencyReferenceLocations, Tuple2<SvType, Iterable<ChimericAlignment>>>
//    inferInvDupRange(final Tuple2<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> noveltyAndEvidence,
//                     final ReferenceMultiSource reference) {
//        return null;
//    }
}
