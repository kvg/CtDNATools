package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;

@DefaultSerializer(EvidenceTargetLink.Serializer.class)
class EvidenceTargetLink {
    private static final SVInterval.Serializer intervalSerializer = new SVInterval.Serializer();

    final SVInterval source;
    final boolean sourceForwardStrand;
    final SVInterval target;
    final boolean targetForwardStrand;
    final int directedWeight;
    final int undirectedWeight;

    public EvidenceTargetLink(final SVInterval source, final boolean sourceForwardStrand, final SVInterval target, final boolean targetForwardStrand, final int directedWeight, final int undirectedWeight) {
        this.source = source;
        this.sourceForwardStrand = sourceForwardStrand;
        this.target = target;
        this.targetForwardStrand = targetForwardStrand;
        this.directedWeight = directedWeight;
        this.undirectedWeight = undirectedWeight;
    }

    public EvidenceTargetLink(final Kryo kryo, final Input input) {
        this.source = intervalSerializer.read(kryo, input, SVInterval.class);
        this.sourceForwardStrand = input.readBoolean();
        this.target = intervalSerializer.read(kryo, input, SVInterval.class);
        this.targetForwardStrand = input.readBoolean();

        this.directedWeight = input.readInt();
        this.undirectedWeight = input.readInt();
    }

    protected void serialize(final Kryo kryo, final Output output ) {
        intervalSerializer.write(kryo, output, source);
        output.writeBoolean(sourceForwardStrand);
        intervalSerializer.write(kryo, output, target);
        output.writeBoolean(targetForwardStrand);

        output.writeInt(directedWeight);
        output.writeInt(undirectedWeight);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<EvidenceTargetLink> {
        @Override
        public void write( final Kryo kryo, final Output output, final EvidenceTargetLink evidence ) {
            evidence.serialize(kryo, output);
        }

        @Override
        public EvidenceTargetLink read(final Kryo kryo, final Input input, final Class<EvidenceTargetLink> klass ) {
            return new EvidenceTargetLink(kryo, input);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final EvidenceTargetLink link = (EvidenceTargetLink) o;

        if (sourceForwardStrand != link.sourceForwardStrand) return false;
        if (targetForwardStrand != link.targetForwardStrand) return false;
        if (directedWeight != link.directedWeight) return false;
        if (undirectedWeight != link.undirectedWeight) return false;
        if (source != null ? !source.equals(link.source) : link.source != null) return false;
        return target != null ? target.equals(link.target) : link.target == null;
    }

    @Override
    public int hashCode() {
        int result = source != null ? source.hashCode() : 0;
        result = 31 * result + (sourceForwardStrand ? 1 : 0);
        result = 31 * result + (target != null ? target.hashCode() : 0);
        result = 31 * result + (targetForwardStrand ? 1 : 0);
        result = 31 * result + directedWeight;
        result = 31 * result + undirectedWeight;
        return result;
    }

    @Override
    public String toString() {
        return "EvidenceTargetLink{" +
                "source=" + source +
                ", sourceForwardStrand=" + sourceForwardStrand +
                ", target=" + target +
                ", targetForwardStrand=" + targetForwardStrand +
                ", directedWeight=" + directedWeight +
                ", undirectedWeight=" + undirectedWeight +
                '}';
    }
}
