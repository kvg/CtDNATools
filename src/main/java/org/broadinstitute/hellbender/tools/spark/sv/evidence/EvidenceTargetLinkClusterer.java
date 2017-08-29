package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervals;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

public class EvidenceTargetLinkClusterer {

    private final ReadMetadata readMetadata;
    private final int minEvidenceMapq;

    public EvidenceTargetLinkClusterer(final ReadMetadata readMetadata, final int minEvidenceMapq) {
        this.readMetadata = readMetadata;
        this.minEvidenceMapq = minEvidenceMapq;
    }

    public Iterator<EvidenceTargetLink> cluster(final Iterator<BreakpointEvidence> breakpointEvidenceIterator) throws Exception {
        final List<EvidenceTargetLink> links = new ArrayList<>();
        final SVIntervalTree<EvidenceTargetLink> currentIntervalsWithTargets = new SVIntervalTree<>();
        while (breakpointEvidenceIterator.hasNext()) {
            final BreakpointEvidence nextEvidence = breakpointEvidenceIterator.next();
            if (nextEvidence.hasDistalTargets(readMetadata, minEvidenceMapq)) {
                Utils.validate(nextEvidence instanceof BreakpointEvidence.SplitRead || nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence,
                        "Unknown evidence type with distal target: " + nextEvidence);
                EvidenceTargetLink updatedLink = null;
                for (final Iterator<SVIntervalTree.Entry<EvidenceTargetLink>> it = currentIntervalsWithTargets.overlappers(nextEvidence.getLocation()); it.hasNext(); ) {
                    final SVIntervalTree.Entry<EvidenceTargetLink> sourceIntervalEntry = it.next();
                    final EvidenceTargetLink oldLink = sourceIntervalEntry.getValue();
                    // todo: what to do if there are more than one distal targets -- for now just taking the first one
                    // this would only be an issue with split reads with more than one SA mapping
                    if (nextEvidence.hasDistalTargets(readMetadata, minEvidenceMapq) &&
                            strandsMatch(nextEvidence.isForwardStrand(), sourceIntervalEntry.getValue().sourceForwardStrand)
                            && (nextEvidence.getDistalTargets(readMetadata, minEvidenceMapq).get(0).overlaps(oldLink.target) &&
                            strandsMatch(nextEvidence.getDistalTargetStrands(readMetadata, minEvidenceMapq).get(0), oldLink.targetForwardStrand))) {
                        // if it does, intersect the source and target intervals to refine the link
                        it.remove();
                        final SVInterval newSource = sourceIntervalEntry.getInterval().intersect(nextEvidence.getLocation());
                        final SVInterval newTarget = oldLink.target.intersect(nextEvidence.getDistalTargets(readMetadata, minEvidenceMapq).get(0));
                        int newSplitReadCount = nextEvidence instanceof BreakpointEvidence.SplitRead
                                ? oldLink.splitReads + 1 : oldLink.splitReads;
                        int newReadPairCount = nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence
                                ? oldLink.readPairs + 1 : oldLink.readPairs;
                        updatedLink = new EvidenceTargetLink(newSource,
                                oldLink.sourceForwardStrand,
                                newTarget,
                                oldLink.targetForwardStrand,
                                newSplitReadCount,
                                newReadPairCount);
                            break;
                    }
                }
                if (updatedLink == null) {
                    updatedLink = new EvidenceTargetLink(
                            nextEvidence.getLocation(),
                            nextEvidence.isForwardStrand(),
                            nextEvidence.getDistalTargets(readMetadata, minEvidenceMapq).get(0),
                            nextEvidence.getDistalTargetStrands(readMetadata, minEvidenceMapq).get(0),
                            nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence
                                    ? 0 : 1,
                            nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence
                                    ? 1 : 0);
                }
                currentIntervalsWithTargets.put(updatedLink.source, updatedLink);
            }
        }

        currentIntervalsWithTargets.forEach(entry -> links.add(entry.getValue()));

        return links.iterator();
    }

    private static boolean strandsMatch(final Boolean forwardStrand1, final Boolean forwardStrand2) {
        return forwardStrand1 != null && forwardStrand2 != null && forwardStrand1.equals(forwardStrand2);
    }

    /**
     * Combines links that agree on interval-pair and orientation but have source and target switched. For
     * example, if link1 has intervals s1 and t1, and link2 has s2 and t2, and s1 overlaps t1 and s2 overlaps t2,
     * the links will be combined as long as strands agree. Returned links have the intersection of two paired intervals as
     * source and target, with the lower-coordinate interval appearing as source.
     */
    public static List<EvidenceTargetLink> deduplicateTargetLinks(final List<EvidenceTargetLink> evidenceTargetLinks) {

        PairedStrandedIntervalTree<EvidenceTargetLink> pairedStrandedIntervalTree = new PairedStrandedIntervalTree<>();

        evidenceTargetLinks.stream().filter(link -> link.source.compareTo(link.target) < 0).forEach(link -> pairedStrandedIntervalTree.put(link.getPairedStrandedIntervals(), link));


        evidenceTargetLinks.stream().filter(link -> link.source.compareTo(link.target) >= 0).forEach(link -> {

            PairedStrandedIntervals reversedPair =
                    new PairedStrandedIntervals(link.target, link.targetForwardStrand, link.source, link.sourceForwardStrand);

            Iterator<Tuple2<PairedStrandedIntervals, EvidenceTargetLink>> psiOverlappers =
                    pairedStrandedIntervalTree.overlappers(reversedPair);

            EvidenceTargetLink newLink = null;
            while (psiOverlappers.hasNext()) {
                EvidenceTargetLink existingLink = psiOverlappers.next()._2();
                newLink = new EvidenceTargetLink(link.target.intersect(existingLink.source), link.targetForwardStrand,
                        link.source.intersect(existingLink.target), link.sourceForwardStrand,
                        link.splitReads + existingLink.splitReads, Math.max(link.readPairs, existingLink.readPairs));
                psiOverlappers.remove();
                break;
            }

            if (newLink != null) {
                pairedStrandedIntervalTree.put(newLink.getPairedStrandedIntervals(), newLink);
            } else {
                pairedStrandedIntervalTree.put(reversedPair, new EvidenceTargetLink(reversedPair.getLeft(), reversedPair.getLeftStrand(), reversedPair.getRight(), reversedPair.getRightStrand(), link.splitReads, link.readPairs));
            }

        });

        return Utils.stream(pairedStrandedIntervalTree.iterator()).map(Tuple2::_2).collect(Collectors.toList());
    }

}
