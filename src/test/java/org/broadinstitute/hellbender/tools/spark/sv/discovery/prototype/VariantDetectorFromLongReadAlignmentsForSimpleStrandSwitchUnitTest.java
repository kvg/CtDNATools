package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.IOException;

public class VariantDetectorFromLongReadAlignmentsForSimpleStrandSwitchUnitTest extends BaseTest {

    @DataProvider(name = "forComputeNewRefSpanAndCigar")
    private Object[][] createTestData() throws IOException {

        final Object[][] data = new Object[4][];

        AlignmentInterval interval = new AlignmentInterval(new SimpleInterval("chr1", 175417007, 175417074),
                14, 81, TextCigarCodec.decode("13H68M394H"),
                true, 60, 0, 68, false, false);
        data[0] = new Object[]{interval, 31, true, new SimpleInterval("chr1", 175417007, 175417043), TextCigarCodec.decode("13H37M31S394H")};

        interval = new AlignmentInterval(new SimpleInterval("chr2", 122612655, 122612751),
                9, 105, TextCigarCodec.decode("8H97M138H"),
                false, 60, 0, 97, false, false);
        data[1] = new Object[]{interval, 4, true, new SimpleInterval("chr2", 122612659, 122612751), TextCigarCodec.decode("8H93M4S138H")};

        interval = new AlignmentInterval(new SimpleInterval("chr6", 66782514, 66782679),
                32, 197, TextCigarCodec.decode("31S166M"),
                false, 60, 3, 151, false, false);
        data[2] = new Object[]{interval, 4, false, new SimpleInterval("chr6", 66782514, 66782675), TextCigarCodec.decode("35S162M")};

        interval = new AlignmentInterval(new SimpleInterval("chr2", 91421528, 91421734),
                271, 477, TextCigarCodec.decode("270H207M"),
                true, 40, 12, 147, false, false);
        data[3] = new Object[]{interval, 32, false, new SimpleInterval("chr2", 91421560, 91421734), TextCigarCodec.decode("270H32S175M")};

        return data;
    }

    @Test(dataProvider = "forComputeNewRefSpanAndCigar", groups = "sv")
    public void testComputeNewRefSpanAndCigar(final AlignmentInterval interval, final int clipLength, final boolean clipFrom3PrimeEnd,
                                              final SimpleInterval expectedRefSpan, final Cigar expectedCigar) {

        final Tuple2<SimpleInterval, Cigar> x = ForSimpleStrandSwitch.computeNewRefSpanAndCigar(interval, clipLength, clipFrom3PrimeEnd);
        Assert.assertEquals(x._1, expectedRefSpan);
        Assert.assertEquals(x._2, expectedCigar);
    }
}
