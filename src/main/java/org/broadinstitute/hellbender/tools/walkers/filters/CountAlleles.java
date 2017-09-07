package org.broadinstitute.hellbender.tools.walkers.filters;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * Created by kiran on 07/09/2017.
 */
@CommandLineProgramProperties(
        summary = "Count alleles",
        oneLineSummary = "Count alleles",
        programGroup = ReadProgramGroup.class
)
@DocumentedFeature
public class CountAlleles extends LocusWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File OUTPUT_FILE = null;

    @Argument(fullName="variants", shortName="V", doc="Set of variants", optional=true)
    private FeatureInput<VariantContext> variants;

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
        if ( featureContext.hasBackingDataSource() ) {
            for (VariantContext vc : featureContext.getValues(variants)) {
                if (vc.getFilters().isEmpty()) {
                    int[] counts = alignmentContext.getBasePileup().getBaseCounts();

                    int altCount = 0;
                    String alt = vc.getAltAlleleWithHighestAlleleCount().getBaseString();
                    if (alt.equals("A")) { altCount = counts[0]; }
                    if (alt.equals("C")) { altCount = counts[1]; }
                    if (alt.equals("G")) { altCount = counts[2]; }
                    if (alt.equals("T")) { altCount = counts[3]; }

                    int totalCount = counts[0] + counts[1] + counts[2] + counts[3];

                    outputStream.println(referenceContext.getInterval().getContig() + " " + referenceContext.getInterval().getStart() + " " + vc.getReference() + " " + vc.getAltAlleleWithHighestAlleleCount() + " " + altCount + " " + totalCount);
                }
            }
        }
    }

    @Override
    public void closeTool() {
        outputStream.close();
    }
}
