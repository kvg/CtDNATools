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
public class CountAlleles extends VariantWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File outputFile = null;

    @Argument(fullName="auxiliaryVariants", shortName="av", doc="Auxiliary set of variants", optional=true)
    private FeatureInput<VariantContext> auxiliaryVariants;

    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        try {
            outputStream = outputFile != null ? new PrintStream(outputFile) : System.out;
        }
        catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(outputFile, e);
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (variant.getFilters().isEmpty()) {
            Allele refAllele = variant.getReference();
            Allele altAllele = variant.getAltAlleleWithHighestAlleleCount();

            int altCount = 0;
            int totalCount = 0;

            int start = variant.getStart();

            for (GATKRead read : readsContext) {
                int offsetStart = start - read.getAssignedStart();
                int offsetEnd = offsetStart + altAllele.length();

                String seenAllele = read.getBasesString().substring(offsetStart, offsetEnd);

                outputStream.println(refAllele + " " + altAllele + " " + seenAllele + " " + variant.getContig() + " " + variant.getStart() + " " + read.getAssignedStart() + " " + offsetStart + " " + read.getBasesString());

                if (altAllele.getDisplayString().equals(seenAllele)) {
                    altCount++;
                }

                totalCount++;
            }

            outputStream.println(altCount + " " + totalCount);
        }

        /*
        outputStream.println("Current variant: " + variant);

        if ( referenceContext.hasBackingDataSource() ) {
            printReferenceBases(referenceContext);
        }

        if ( readsContext.hasBackingDataSource() ) {
            printReads(readsContext);
        }

        if ( featureContext.hasBackingDataSource() ) {
            printVariants(featureContext);
        }
        */
    }

    private void printReferenceBases( final ReferenceContext refContext ) {
        outputStream.printf("\tOverlapping reference bases: %s\n\n", new String(refContext.getBases()));
    }

    private void printReads( final ReadsContext readsContext ) {
        for ( final GATKRead read : readsContext ) {
            outputStream.printf("\tOverlapping read at %s:%d-%d\n", read.getContig(), read.getStart(), read.getEnd());
        }
        outputStream.println();
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
        if ( outputStream != null ) {
            outputStream.close();
        }
    }
}
