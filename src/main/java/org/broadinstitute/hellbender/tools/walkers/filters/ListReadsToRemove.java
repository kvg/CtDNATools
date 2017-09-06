package org.broadinstitute.hellbender.tools.walkers.filters;

import ngs.Pileup;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by kiran on 06/09/2017.
 */
@CommandLineProgramProperties(
        summary = "Filter out reads that are the sole evidence for mismatches",
        oneLineSummary = "Filter out reads that are the sole evidence for mismatches",
        programGroup = ReadProgramGroup.class
)
@DocumentedFeature
public class ListReadsToRemove extends LocusWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File OUTPUT_FILE = null;

    private Set<GATKRead> toRemove = new HashSet<>();
    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        try {
            outputStream = OUTPUT_FILE != null ? new PrintStream(OUTPUT_FILE) : System.out;
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(OUTPUT_FILE, e);
        }
    }

    @Override
    public void apply(AlignmentContext alignmentContext, final ReferenceContext referenceContext, FeatureContext featureContext) {
        ReadPileup rp = alignmentContext.getBasePileup().makeFilteredPileup(pe -> pe.getBase() != referenceContext.getBase());

        if (rp.getReads().size() == 1) {
            toRemove.add(rp.getReads().get(0));
        }
    }

    @Override
    public void closeTool() {
        for (GATKRead gr : toRemove) {
            outputStream.println(gr.getName());
        }

        outputStream.close();
    }
}
