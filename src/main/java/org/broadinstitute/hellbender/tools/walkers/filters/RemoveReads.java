package org.broadinstitute.hellbender.tools.walkers.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;

import java.io.*;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by kiran on 06/09/2017.
 */
@CommandLineProgramProperties(
        summary = "Filter out specified reads",
        oneLineSummary = "Filter out specified reads",
        programGroup = ReadProgramGroup.class
)
@DocumentedFeature
public class RemoveReads extends ReadWalker {
    @Argument(fullName="readList", shortName="rl", doc="Read list")
    public File READ_LIST;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file")
    public File OUTPUT;

    private SAMFileGATKReadWriter outputWriter;
    private Set<String> toRemove = new HashSet<>();

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(OUTPUT, true);

        try {
            FileInputStream fis = new FileInputStream(READ_LIST);
            BufferedReader reader = new BufferedReader(new InputStreamReader(fis));

            String line = reader.readLine();
            while (line != null) {
                line = reader.readLine();

                toRemove.add(line);
            }
        } catch (FileNotFoundException e) {
            throw new GATKException("File not found", e);
        } catch (IOException e) {
            throw new GATKException("IO exception", e);
        }
    }

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        if (!toRemove.contains(read.getName())) {
            outputWriter.addRead(read);
        }
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }
}
