package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.spark.api.java.function.FlatMapFunction;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class EvidenceTargetLinkClusterer implements FlatMapFunction<Iterator<BreakpointEvidence>, EvidenceTargetLink> {

    private final ReadMetadata readMetadata;
    private final int totalNumIntervals;

    public EvidenceTargetLinkClusterer(final ReadMetadata readMetadata, final int totalNumIntervals) {
        this.readMetadata = readMetadata;
        this.totalNumIntervals = totalNumIntervals;
    }

    @Override
    public Iterator<EvidenceTargetLink> call(final Iterator<BreakpointEvidence> breakpointEvidenceIterator) throws Exception {
        final List<EvidenceTargetLink> links = new ArrayList<>(totalNumIntervals / readMetadata.getNPartitions());
        final SVIntervalTree<EvidenceTargetLink> currentIntervalsWithTargets = new SVIntervalTree<>();
        while (breakpointEvidenceIterator.hasNext()) {
            final BreakpointEvidence nextEvidence = breakpointEvidenceIterator.next();
            if (nextEvidence.hasDistalTargets()) {
                EvidenceTargetLink updatedLink = null;
                for (final Iterator<SVIntervalTree.Entry<EvidenceTargetLink>> it = currentIntervalsWithTargets.iterator(); it.hasNext(); ) {
                    final SVIntervalTree.Entry<EvidenceTargetLink> sourceIntervalEntry = it.next();
                    final EvidenceTargetLink oldLink = sourceIntervalEntry.getValue();
                    // todo: what to do if there are more than one distal targets
                    if (nextEvidence.hasDistalTargets() && nextEvidence.getLocation().overlaps(sourceIntervalEntry.getInterval())
                            && strandsMatch(nextEvidence.isForwardStrand(), sourceIntervalEntry.getValue().sourceForwardStrand)
                        && (nextEvidence.getDistalTargets(readMetadata).get(0).overlaps(oldLink.target) &&
                                strandsMatch(nextEvidence.getDistalTargetStrands().get(0), oldLink.targetForwardStrand))) {
                            // if it does, intersect the source and target intervals to refine the link
                            it.remove();
                            final SVInterval newSource = sourceIntervalEntry.getInterval().intersect(nextEvidence.getLocation());
                            final SVInterval newTarget = oldLink.target.intersect(nextEvidence.getDistalTargets(readMetadata).get(0));
                            updatedLink = new EvidenceTargetLink(newSource,
                                    oldLink.sourceForwardStrand,
                                    newTarget,
                                    oldLink.targetForwardStrand,
                                    nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence
                                            ? oldLink.directedWeight : oldLink.directedWeight + 1,
                                    nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence
                                            ? oldLink.undirectedWeight + 1 : oldLink.undirectedWeight);
                            break;
                    }
                }
                if (updatedLink == null) {
                    updatedLink = new EvidenceTargetLink(
                            nextEvidence.getLocation(),
                            nextEvidence.isForwardStrand(),
                            nextEvidence.getDistalTargets(readMetadata).get(0),
                            nextEvidence.getDistalTargetStrands().get(0),
                            nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence
                                    ? 0 : 1,
                            nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence
                                    ? 1 : 0);

                }
                currentIntervalsWithTargets.put(updatedLink.source, updatedLink);
            }
        }
        for (Iterator<SVIntervalTree.Entry<EvidenceTargetLink>> it = currentIntervalsWithTargets.iterator(); it.hasNext(); ) {
            links.add(it.next().getValue());
        }
        return links.iterator();
    }

    private static boolean strandsMatch(final Boolean forwardStrand1, final Boolean forwardStrand2) {
        return forwardStrand1 != null && forwardStrand2 != null && forwardStrand1.equals(forwardStrand2);
    }

}
