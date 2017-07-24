package org.broadinstitute.hellbender.tools.spark.sv;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class EvidenceTargetLinkClusterer {

    private final ReadMetadata readMetadata;

    public EvidenceTargetLinkClusterer(final ReadMetadata readMetadata) {
        this.readMetadata = readMetadata;
    }

    public Iterator<EvidenceTargetLink> cluster(final Iterator<BreakpointEvidence> breakpointEvidenceIterator) throws Exception {
        final List<EvidenceTargetLink> links = new ArrayList<>();
        final SVIntervalTree<EvidenceTargetLink> currentIntervalsWithTargets = new SVIntervalTree<>();
        while (breakpointEvidenceIterator.hasNext()) {
            final BreakpointEvidence nextEvidence = breakpointEvidenceIterator.next();
            if (nextEvidence.hasDistalTargets(readMetadata)) {
                final SVInterval target = nextEvidence.getDistalTargets(readMetadata).get(0);
                System.err.println(toBedpeString(nextEvidence, nextEvidence.getLocation(), target, ((BreakpointEvidence.ReadEvidence) nextEvidence).getTemplateName() +
                        ((BreakpointEvidence.ReadEvidence) nextEvidence).getFragmentOrdinal() + (nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence ? "DRP" : "SR"), 1));
                EvidenceTargetLink updatedLink = null;
                for (final Iterator<SVIntervalTree.Entry<EvidenceTargetLink>> it = currentIntervalsWithTargets.overlappers(nextEvidence.getLocation()); it.hasNext(); ) {
                    final SVIntervalTree.Entry<EvidenceTargetLink> sourceIntervalEntry = it.next();
                    final EvidenceTargetLink oldLink = sourceIntervalEntry.getValue();
                    // todo: what to do if there are more than one distal targets
                    if (nextEvidence.hasDistalTargets(readMetadata) &&
                            strandsMatch(nextEvidence.isForwardStrand(), sourceIntervalEntry.getValue().sourceForwardStrand)
                        && (nextEvidence.getDistalTargets(readMetadata).get(0).overlaps(oldLink.target) &&
                                strandsMatch(nextEvidence.getDistalTargetStrands(readMetadata).get(0), oldLink.targetForwardStrand))) {
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
                            System.err.println("updating to: " + toBedpeString(nextEvidence, newSource, newTarget, ((BreakpointEvidence.ReadEvidence) nextEvidence).getTemplateName() +
                                    ((BreakpointEvidence.ReadEvidence) nextEvidence).getFragmentOrdinal() + "_" + updatedLink.directedWeight + "_" + updatedLink.undirectedWeight, 1));
                            break;
                    }
                }
                if (updatedLink == null) {
                    updatedLink = new EvidenceTargetLink(
                            nextEvidence.getLocation(),
                            nextEvidence.isForwardStrand(),
                            nextEvidence.getDistalTargets(readMetadata).get(0),
                            nextEvidence.getDistalTargetStrands(readMetadata).get(0),
                            nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence
                                    ? 0 : 1,
                            nextEvidence instanceof BreakpointEvidence.DiscordantReadPairEvidence
                                    ? 1 : 0);
                    System.err.println("creating new: " + toBedpeString(nextEvidence, nextEvidence.getLocation(), nextEvidence.getDistalTargets(readMetadata).get(0), ((BreakpointEvidence.ReadEvidence) nextEvidence).getTemplateName() +
                            ((BreakpointEvidence.ReadEvidence) nextEvidence).getFragmentOrdinal() + "_" + updatedLink.directedWeight + "_" + updatedLink.undirectedWeight, 1));
                }
                currentIntervalsWithTargets.put(updatedLink.source, updatedLink);
            }
        }
        for (Iterator<SVIntervalTree.Entry<EvidenceTargetLink>> it = currentIntervalsWithTargets.iterator(); it.hasNext(); ) {
            links.add(it.next().getValue());
        }
        return links.iterator();
    }

    public String toBedpeString(final BreakpointEvidence nextEvidence, final SVInterval source, final SVInterval target, final String id, final int score) {
        return "21" + "\t" + (source.getStart() - 1) + "\t" + source.getEnd() +
                "\t" + "21" + "\t" + (target.getStart() - 1) + "\t" + target.getEnd() + "\t"  +
                id + "\t" +
                score + "\t" + (nextEvidence.isForwardStrand() ? "+" : "-") + "\t" + (nextEvidence.getDistalTargetStrands(readMetadata).get(0) ? "+" : "-");
    }

    private static boolean strandsMatch(final Boolean forwardStrand1, final Boolean forwardStrand2) {
        return forwardStrand1 != null && forwardStrand2 != null && forwardStrand1.equals(forwardStrand2);
    }

}
