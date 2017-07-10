package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.StringUtil;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class TextMDCodec {

    /**
     * Regexp for MD string. Copied from htsjdk.samtools.util.SequenceUtil
     * //todo: if this parser is useful, make a pull request to add it to htsjdk
     *
     * \G = end of previous match.
     * (?:[0-9]+) non-capturing (why non-capturing?) group of digits.  For this number of bases read matches reference.
     *  - or -
     * Single reference base for case in which reference differs from read.
     *  - or -
     * ^one or more reference bases that are deleted in read.
     *
     */
    static final Pattern mdPat = Pattern.compile("\\G(?:([0-9]+)|([ACTGNactgn])|(\\^[ACTGNactgn]+))");

    static List<MDElement> parseMDString(final String mdString) {
        List<MDElement> results = new ArrayList<>();
        final Matcher match = mdPat.matcher(mdString);
        while (match.find()) {
            String mg;
            if (((mg = match.group(1)) != null) && (!mg.isEmpty())) {
                // It's a number , meaning a series of matches
                final int num = Integer.parseInt(mg);
                results.add(new MatchMDElement(num));
            } else if (((mg = match.group(2)) != null) && (!mg.isEmpty())) {
                results.add(new MismatchMDElement());
            } else if (((mg = match.group(3)) != null) && (!mg.isEmpty())) {
                // remove the carat and the deletion length it what's left
                results.add(new DeletionMDElement(mg.length() - 1));
            }
        }
        return results;
    }

    abstract static class MDElement {
        public abstract int getLength();
    }

    static class MatchMDElement extends MDElement {

        final int length;

        public MatchMDElement(final int length) {
            this.length = length;
        }

        @Override
        public int getLength() {
            return length;
        }
    }

    static class MismatchMDElement extends MDElement {

        @Override
        public int getLength() {
            return 1;
        }
    }

    static class DeletionMDElement extends MDElement {

        final int length;

        public DeletionMDElement(final int length) {
            this.length = length;
        }

        @Override
        public int getLength() {
            return length;
        }
    }
}
