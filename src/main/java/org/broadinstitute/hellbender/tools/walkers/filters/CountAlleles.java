package org.broadinstitute.hellbender.tools.walkers.filters;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;

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

    @Argument(fullName="auxiliaryVariants", shortName="av", doc="Auxiliary set of variants", optional=true)
    private FeatureInput<VariantContext> auxiliaryVariants;

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
            printVariants(featureContext);
        }

        int[] counts = alignmentContext.getBasePileup().getBaseCounts();

        outputStream.println(referenceContext.getInterval().getContig() + " " + referenceContext.getInterval().getStart() + " " + counts[0] + " " + counts[1] + " " + counts[2] + " " + counts[3]);
    }

    private void printVariants( final FeatureContext featureContext ) {
        for ( final VariantContext variant : featureContext.getValues(auxiliaryVariants) ) {
            outputStream.printf("\tOverlapping variant at %s:%d-%d. Ref: %s Alt(s): %s\n",
                    variant.getContig(), variant.getStart(), variant.getEnd(), variant.getReference(), variant.getAlternateAlleles());
        }
        outputStream.println();
    }

    @Override
    public void closeTool() {
        outputStream.close();
    }
}
