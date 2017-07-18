package org.broadinstitute.hellbender.tools.spark.sv;

import javafx.util.Pair;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class PairedStrandedIntervalTree {

    private SVIntervalTree<Tuple2<Boolean, SVIntervalTree<Boolean>>> leftEnds = new SVIntervalTree<>();

    public boolean put(PairedStrandedIntervals pair) {
        if (contains(pair)) return false;

        final SVIntervalTree.Entry<Tuple2<Boolean, SVIntervalTree<Boolean>>> leftEntry = leftEnds.find(pair.getLeft());
        if (leftEntry != null) {
            leftEntry.getValue()._2().put(pair.right, pair.rightStrand);
        } else {
            final SVIntervalTree<Boolean> rightEnds = new SVIntervalTree<>();
            rightEnds.put(pair.right, pair.getRightStrand());
            leftEnds.put(pair.getLeft(), new Tuple2<>(pair.getLeftStrand(), rightEnds));
        }

        return true;
    }

    public Iterator<PairedStrandedIntervals> overlappers(PairedStrandedIntervals pair) {
        final List<PairedStrandedIntervals> pairedStrandedIntervalsList = new ArrayList<>();
        final Iterator<SVIntervalTree.Entry<Tuple2<Boolean, SVIntervalTree<Boolean>>>> leftOverlappers = leftEnds.overlappers(pair.getLeft());
        while (leftOverlappers.hasNext()) {
            SVIntervalTree.Entry<Tuple2<Boolean, SVIntervalTree<Boolean>>> leftEntry = leftOverlappers.next();
            final Boolean leftEntryStrand = leftEntry.getValue()._1();
            if (leftEntryStrand == pair.getLeftStrand()) {
                final SVIntervalTree<Boolean> rightIntervals = leftEntry.getValue()._2();
                final Iterator<SVIntervalTree.Entry<Boolean>> rightOverlappers = rightIntervals.overlappers(pair.getRight());
                while (rightOverlappers.hasNext()) {
                    SVIntervalTree.Entry<Boolean> rightEntry = rightOverlappers.next();
                    if (rightEntry.getValue() == pair.getRightStrand()) {
                        pairedStrandedIntervalsList.add(new PairedStrandedIntervals(leftEntry.getInterval(), leftEntryStrand, rightEntry.getInterval(), rightEntry.getValue()));
                    }
                }
            }
        }
        return Collections.unmodifiableList(pairedStrandedIntervalsList).iterator();
    }

    public Iterator<PairedStrandedIntervals> iterator() {
        final List<PairedStrandedIntervals> pairedStrandedIntervalsList = new ArrayList<>(leftEnds.size() * 2);
        final Iterator<SVIntervalTree.Entry<Tuple2<Boolean, SVIntervalTree<Boolean>>>> leftEndIterator = leftEnds.iterator();
        while (leftEndIterator.hasNext()) {
            final SVIntervalTree.Entry<Tuple2<Boolean, SVIntervalTree<Boolean>>> leftEntry = leftEndIterator.next();
            final Boolean leftEntryStrand = leftEntry.getValue()._1();
            final Iterator<SVIntervalTree.Entry<Boolean>> rightEndIterator = leftEntry.getValue()._2().iterator();
            while (rightEndIterator.hasNext()) {
                SVIntervalTree.Entry<Boolean> rightEntry = rightEndIterator.next();
                pairedStrandedIntervalsList.add(new PairedStrandedIntervals(leftEntry.getInterval(), leftEntryStrand, rightEntry.getInterval(), rightEntry.getValue()));
            }
        }
        return Collections.unmodifiableList(pairedStrandedIntervalsList).iterator();
    }

    public boolean contains(PairedStrandedIntervals pair) {
        final int leftEndIndex = leftEnds.getIndex(pair.getLeft());
        if (leftEndIndex == -1) return false;
        final SVIntervalTree.Entry<Tuple2<Boolean, SVIntervalTree<Boolean>>> leftEndEntry = leftEnds.findByIndex(leftEndIndex);
        final Tuple2<Boolean, SVIntervalTree<Boolean>> storedValue = leftEndEntry.getValue();

        if (pair.leftStrand != storedValue._1()) return false;

        final SVIntervalTree<Boolean> rightEnds = storedValue._2();
        final int rightIndex = rightEnds.getIndex(pair.getRight());
        return rightIndex != -1 && (pair.rightStrand == rightEnds.findByIndex(rightIndex).getValue());
    }

    public static final class PairedStrandedIntervals {
        SVInterval left;
        boolean leftStrand;
        SVInterval right;
        boolean rightStrand;

        public PairedStrandedIntervals(final SVInterval left, final boolean leftStrand, final SVInterval right, final boolean rightStrand) {
            this.left = left;
            this.leftStrand = leftStrand;
            this.right = right;
            this.rightStrand = rightStrand;
        }

        public SVInterval getLeft() {
            return left;
        }

        public void setLeft(final SVInterval left) {
            this.left = left;
        }

        public boolean getLeftStrand() {
            return leftStrand;
        }

        public void setLeftStrand(final boolean leftStrand) {
            this.leftStrand = leftStrand;
        }

        public SVInterval getRight() {
            return right;
        }

        public void setRight(final SVInterval right) {
            this.right = right;
        }

        public boolean getRightStrand() {
            return rightStrand;
        }

        public void setRightStrand(final boolean rightStrand) {
            this.rightStrand = rightStrand;
        }
    }
}
