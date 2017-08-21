package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.Iterator;
import java.util.function.BiPredicate;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection;

@DefaultSerializer(SVReadFilter.Serializer.class)
public class SVReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;

    private final int minEvidenceMapQ;
    private final int minEvidenceMatchLength;
    private final int allowedShortFragmentOverhang;

    public SVReadFilter( final FindBreakpointEvidenceSparkArgumentCollection params ) {
        minEvidenceMapQ = params.minEvidenceMapQ;
        minEvidenceMatchLength = params.minEvidenceMatchLength;
        allowedShortFragmentOverhang = params.allowedShortFragmentOverhang;
    }

    public SVReadFilter(final Kryo kryo, final Input input) {
        this.minEvidenceMapQ = input.readInt();
        this.minEvidenceMatchLength = input.readInt();
        this.allowedShortFragmentOverhang = input.readInt();
    }

    public boolean notJunk( final GATKRead read ) {
        return !read.isDuplicate() && !read.failsVendorQualityCheck();
    }

    public boolean isPrimaryLine( final GATKRead read ) {
        return !read.isSecondaryAlignment() && !read.isSupplementaryAlignment();
    }

    public boolean isMapped( final GATKRead read ) {
        return notJunk(read) && !read.isUnmapped();
    }

    public boolean isMappedPrimary( final GATKRead read ) {
        return isMapped(read) && isPrimaryLine(read);
    }

    public boolean isEvidence( final GATKRead read ) {
        return isMapped(read) && read.getMappingQuality() >= minEvidenceMapQ &&
                SVUtils.matchLen(read.getCigar()) >= minEvidenceMatchLength && ! read.isSecondaryAlignment();
    }

    public boolean isTemplateLenTestable( final GATKRead read ) {
        return isEvidence(read) && isPrimaryLine(read) &&
                !read.mateIsUnmapped() &&
                !read.isReverseStrand() &&
                read.mateIsReverseStrand() &&
                read.getContig().equals(read.getMateContig()) &&
                read.getStart() - allowedShortFragmentOverhang <= read.getMateStart();
    }

    public Iterator<GATKRead> applyFilter( final Iterator<GATKRead> readItr, final BiPredicate<SVReadFilter, GATKRead> predicate ) {
        return new SVUtils.IteratorFilter<>(readItr, read -> predicate.test(this, read));
    }

    public static class Serializer extends com.esotericsoftware.kryo.Serializer<SVReadFilter> {
        @Override
        public void write(final Kryo kryo, final Output output, final SVReadFilter svReadFilter ) {
            svReadFilter.serialize(kryo, output);
        }

        @Override
        public SVReadFilter read(final Kryo kryo, final Input input, final Class<SVReadFilter> klass ) {
            return new SVReadFilter(kryo, input);
        }

    }

    private void serialize(final Kryo kryo, final Output output) {
        output.writeInt(minEvidenceMapQ);
        output.writeInt(minEvidenceMatchLength);
        output.writeInt(allowedShortFragmentOverhang);
    }

    public int getMinEvidenceMapQ() {
        return minEvidenceMapQ;
    }

    public int getMinEvidenceMatchLength() {
        return minEvidenceMatchLength;
    }

    public int getAllowedShortFragmentOverhang() {
        return allowedShortFragmentOverhang;
    }
}
