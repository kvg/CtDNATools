package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.utils.IntHistogramTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class BreakpointEvidenceTest extends BaseTest {
    private final static LibraryStatistics stats =
            new LibraryStatistics(IntHistogramTest.genLogNormalSample(400, 175, 10000).getCDF(),
                    60000000000L, 600000000L, 3000000000L);

    @Test(groups = "sv")
    void restOfFragmentSizeReverseReadTest() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(1, 1, 10000000, 1);
        final String groupName = header.getReadGroups().get(0).getReadGroupId();
        final int readSize = 151;
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), header, stats, null, 1L, 1L, 1);
        final String templateName = "xyzzy";
        final int readStart = 1010101;
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, templateName, 0, readStart, readSize);
        read.setIsPaired(false);
        read.setIsReverseStrand(true);
        read.setReadGroup(groupName);
        final BreakpointEvidence.ReadEvidence evidence1 = new BreakpointEvidence.ReadEvidence(read, readMetadata);
        final int evidenceWidth = readMetadata.getFragmentLengthStatistics(groupName).getMaxNonOutlierFragmentSize() - readSize;
        // otherwise the test below will break as it is currently structured
        Assert.assertTrue(evidenceWidth % 2 == 0);
        final int uncertainty = evidenceWidth /2;
        final int evidenceLocus = readStart - uncertainty;
        final BreakpointEvidence evidence2 =
                new BreakpointEvidence.ReadEvidence(read, readMetadata, evidenceLocus, uncertainty, ! read.isReverseStrand());
        Assert.assertEquals(evidence1.getLocation(), new SVInterval(0,evidenceLocus-uncertainty,evidenceLocus+uncertainty));
        Assert.assertEquals(evidence1.getLocation().getLength(), 2*uncertainty);
        Assert.assertEquals(evidence1.getTemplateName(), templateName);
        Assert.assertEquals(evidence1.getFragmentOrdinal(), TemplateFragmentOrdinal.UNPAIRED);
        Assert.assertEquals(evidence1.toString(), evidence2.toString());
        read.setIsReverseStrand(false);
        final BreakpointEvidence evidence3 = new BreakpointEvidence.ReadEvidence(read, readMetadata);
        final BreakpointEvidence evidence4 =
                new BreakpointEvidence.ReadEvidence(read, readMetadata, readStart+readSize+uncertainty, uncertainty, ! read.isReverseStrand());
        Assert.assertEquals(evidence3.toString(), evidence4.toString());
    }

    @Test(groups = "sv")
    void restOfFragmentSizeForwardReadTest() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeaderWithGroups(1, 1, 10000000, 1);
        final String groupName = header.getReadGroups().get(0).getReadGroupId();
        final int readSize = 151;
        final ReadMetadata readMetadata = new ReadMetadata(Collections.emptySet(), header, stats, null, 1L, 1L, 1);
        final String templateName = "xyzzy";
        final int readStart = 1010101;
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, templateName, 0, readStart, readSize);
        read.setIsPaired(false);
        read.setIsReverseStrand(false);
        read.setReadGroup(groupName);
        final BreakpointEvidence.ReadEvidence evidence1 = new BreakpointEvidence.ReadEvidence(read, readMetadata);
        final int evidenceWidth = readMetadata.getFragmentLengthStatistics(groupName).getMaxNonOutlierFragmentSize() - readSize;
        // otherwise the test below will break as it is currently structured
        Assert.assertTrue(evidenceWidth % 2 == 0);
        final int uncertainty = evidenceWidth /2;
        final int evidenceLocus = read.getEnd() + 1 + uncertainty;
        final BreakpointEvidence evidence2 =
                new BreakpointEvidence.ReadEvidence(read, readMetadata, evidenceLocus, uncertainty, true);
        Assert.assertEquals(evidence1.getLocation(), new SVInterval(0,evidenceLocus-uncertainty,evidenceLocus+uncertainty));
        Assert.assertEquals(evidence1.getLocation().getLength(), 2*uncertainty);
        Assert.assertEquals(evidence1.getTemplateName(), templateName);
        Assert.assertEquals(evidence1.getFragmentOrdinal(), TemplateFragmentOrdinal.UNPAIRED);
        Assert.assertEquals(evidence1.toString(), evidence2.toString());
        read.setIsReverseStrand(false);
        final BreakpointEvidence evidence3 = new BreakpointEvidence.ReadEvidence(read, readMetadata);
        final BreakpointEvidence evidence4 =
                new BreakpointEvidence.ReadEvidence(read, readMetadata, readStart+readSize+uncertainty, uncertainty, ! read.isReverseStrand());
        Assert.assertEquals(evidence3.toString(), evidence4.toString());
    }

    @Test(groups = "sv")
    void serializationTest() {
        final List<BreakpointEvidence> evidenceList = new ArrayList<>(7);
        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader();
        final ReadMetadata metadata = new ReadMetadata(Collections.emptySet(), samHeader, stats, null, 2L, 2L, 1);
        final List<GATKRead> readPair = ArtificialReadUtils.createPair(samHeader, "firstReadPair", 101, 1010, 1382, false, false);
        final GATKRead read = readPair.get(0);
        evidenceList.add(new BreakpointEvidence.SplitRead(read, metadata, true));
        evidenceList.add(new BreakpointEvidence.LargeIndel(read, metadata, read.getStart()+50));
        evidenceList.add(new BreakpointEvidence.MateUnmapped(read, metadata));
        evidenceList.add(new BreakpointEvidence.InterContigPair(read, metadata));
        evidenceList.add(new BreakpointEvidence.OutiesPair(read, metadata));
        evidenceList.add(new BreakpointEvidence.SameStrandPair(read, metadata));
        evidenceList.add(new BreakpointEvidence.WeirdTemplateSize(read, metadata));

        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        final Output out = new Output(bos);
        final Kryo kryo = new Kryo();
        kryo.writeClassAndObject(out, evidenceList);
        out.flush();

        final ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
        final Input in = new Input(bis);
        @SuppressWarnings("unchecked")
        final List<BreakpointEvidence> evidenceList2 = (List<BreakpointEvidence>)kryo.readClassAndObject(in);
        Assert.assertEquals(evidenceList.size(), evidenceList2.size());
        for ( int idx = 0; idx != evidenceList.size(); ++idx ) {
            Assert.assertEquals(evidenceList.get(idx).toString(), evidenceList2.get(idx).toString());
        }
    }

    @Test(groups = "sv")
    public void testSplitReadFromPrimaryFirstInPair() {
        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader();
        final ReadMetadata metadata = new ReadMetadata(Collections.emptySet(), samHeader, stats, null, 2L, 2L, 1);
        final List<GATKRead> readPair = ArtificialReadUtils.createPair(samHeader, "firstReadPair", 151, 140825480, 140828201, true, false);
        final GATKRead read = readPair.get(0);
        read.setAttribute("SA", "1,140828201,+,82S69M,60,1;");

        final BreakpointEvidence.SplitRead splitRead = new BreakpointEvidence.SplitRead(read, metadata, false);
        Assert.assertTrue(splitRead.isForwardStrand());
        Assert.assertTrue(splitRead.hasDistalTargets(metadata));
        final SVInterval targetInterval = splitRead.getDistalTargets(metadata).get(0);
        Assert.assertEquals(targetInterval, new SVInterval(0, 140828198, 140828274));
        Assert.assertFalse(splitRead.getDistalTargetStrands(metadata).get(0));
    }

    @Test(groups = "sv")
    public void testSplitReadFromPrimarySecondInPair() {
        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader();
        final ReadMetadata metadata = new ReadMetadata(Collections.emptySet(), samHeader, stats, null, 2L, 2L, 1);
        final List<GATKRead> readPair = ArtificialReadUtils.createPair(samHeader, "firstReadPair", 151, 140825480, 140828201, true, false);
        final GATKRead read = readPair.get(1);
        read.setAttribute("SA", "1,140825513,-,20S54M77S,60,0;");

        final BreakpointEvidence.SplitRead splitRead = new BreakpointEvidence.SplitRead(read, metadata, true);
        Assert.assertFalse(splitRead.isForwardStrand());
        Assert.assertTrue(splitRead.hasDistalTargets(metadata));
        final SVInterval targetInterval = splitRead.getDistalTargets(metadata).get(0);
        Assert.assertEquals(targetInterval, new SVInterval(0, 140825510, 140825571));
        Assert.assertTrue(splitRead.getDistalTargetStrands(metadata).get(0));
    }

    @Test(groups = "sv")
    public void testSplitReadFromSupplementaryLeft() {
        final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader();
        final ReadMetadata metadata = new ReadMetadata(Collections.emptySet(), samHeader, stats, null, 2L, 2L, 1);
        final GATKRead read = ArtificialReadUtils.createArtificialRead(samHeader, "saRead", 0, 140828201,
                ArtificialReadUtils.createRandomReadBases(151, false), ArtificialReadUtils.createRandomReadQuals(151), "82S69M");
        read.setIsSupplementaryAlignment(true);
        read.setAttribute("SA", "1,140825480,+,87M64S,60,0;");

        final BreakpointEvidence.SplitRead splitRead = new BreakpointEvidence.SplitRead(read, metadata, true);
        Assert.assertFalse(splitRead.isForwardStrand());
        Assert.assertTrue(splitRead.hasDistalTargets(metadata));
        final SVInterval targetInterval = splitRead.getDistalTargets(metadata).get(0);
        Assert.assertEquals(targetInterval, new SVInterval(0, 140825477, 140825571));
        Assert.assertTrue(splitRead.getDistalTargetStrands(metadata).get(0));
    }

    @Test(groups = "sv")
    public void testSRDistalTargetStrand() {
        Assert.assertFalse(BreakpointEvidence.SplitRead.calculateDistalTargetStrand(
                new BreakpointEvidence.SplitRead.SAMapping("3", 140828201, true, "82S69M", 60, 1), true));

        Assert.assertTrue(BreakpointEvidence.SplitRead.calculateDistalTargetStrand(
                new BreakpointEvidence.SplitRead.SAMapping("3", 140825480, true, "87M64S", 60, 0), false));

        Assert.assertTrue(BreakpointEvidence.SplitRead.calculateDistalTargetStrand(
                new BreakpointEvidence.SplitRead.SAMapping("3", 140825513, false, "20S54M77S", 60, 0), false));

        Assert.assertFalse(BreakpointEvidence.SplitRead.calculateDistalTargetStrand(
                new BreakpointEvidence.SplitRead.SAMapping("3", 140828201, false, "69S82M", 60, 0), true));

        Assert.assertTrue(BreakpointEvidence.SplitRead.calculateDistalTargetStrand(
                new BreakpointEvidence.SplitRead.SAMapping("3", 43593545, false, "81M70S", 60, 2), true));

    }
}
