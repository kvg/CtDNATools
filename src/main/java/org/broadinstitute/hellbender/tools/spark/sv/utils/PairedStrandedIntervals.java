package org.broadinstitute.hellbender.tools.spark.sv.utils;

public final class PairedStrandedIntervals {
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

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final PairedStrandedIntervals that = (PairedStrandedIntervals) o;

        if (leftStrand != that.leftStrand) return false;
        if (rightStrand != that.rightStrand) return false;
        if (left != null ? !left.equals(that.left) : that.left != null) return false;
        return right != null ? right.equals(that.right) : that.right == null;
    }

    @Override
    public int hashCode() {
        int result = left != null ? left.hashCode() : 0;
        result = 31 * result + (leftStrand ? 1 : 0);
        result = 31 * result + (right != null ? right.hashCode() : 0);
        result = 31 * result + (rightStrand ? 1 : 0);
        return result;
    }
}
