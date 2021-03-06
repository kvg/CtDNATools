package org.broadinstitute.hellbender.tools.coveragemodel;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.analysis.interpolation.BicubicInterpolatingFunction;
import org.apache.commons.math3.analysis.interpolation.BicubicInterpolator;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.utils.Utils;
import org.broadinstitute.hellbender.tools.coveragemodel.interfaces.TargetLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.io.InputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

/**
 * This class implements {@link TargetLikelihoodCalculator} according to the GATK Bayesian coverage model.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class CoverageModelCopyRatioEmissionProbabilityCalculator implements
        TargetLikelihoodCalculator<CoverageModelCopyRatioEmissionData>, Serializable {

    private static final long serialVersionUID = -6985799468753075235L;

    /**
     * The following field is for debugging purposes
     */
    private static final boolean CHECK_FOR_NANS = true;

    /**
     * A resource file that provides the range of "mu"s on which we have exact values of
     * the (log) normalization factor of emission probability as a 1D array
     */
    private static final String MU_TABLE_RESOURCE = "mu_table.tsv";

    /**
     * A resource file that provides the range of "psi"s on which we have exact values of
     * the (log) normalization factor of emission probability as a 1D array
     */
    private static final String PSI_TABLE_RESOURCE = "psi_table.tsv";

    /**
     * A resource file that provides exact values of the log normalization factor of
     * the emission probability for (mu, psi) for mu in {@link #MU_TABLE_RESOURCE}
     * and for psi in {@link #PSI_TABLE_RESOURCE}, as a 2D array.
     *
     * These values are calculated via brute-force summation over read count using Mathematica.
     * This table is used for estimating the log normalization factor for a given (mu, psi) pair
     * via bicubic interpolation; see {@link #getLogNormFactorInterpolatingFunction}
     * and {@link #getLogProbabilityLaplaceApproximationNormalizationConstant}.
     */
    private static final String LOG_NORM_TABLE_RESOURCE = "log_norm_table.tsv";

    private static final BicubicInterpolatingFunction logNormFactorSpline = getLogNormFactorInterpolatingFunction();

    /**
     * The following four fields represents min/max values of $\mu$ and $\psi$ in {@link CoverageModelCopyRatioEmissionData}
     * for which the normalization constants of the emission probability distribution are tabulated.
     */
    private static final double MINIMUM_TABULATED_PSI = 0.0;
    private static final double MAXIMUM_TABULATED_PSI = 1.0;
    private static final double MINIMUM_TABULATED_MU = -25.0;
    private static final double MAXIMUM_TABULATED_MU = 10.0;

    private static final double LOG_2PI = FastMath.log(2 * FastMath.PI);

    /**
     * The interpolation scheme used for calculating the log normalization factor of the emission probability
     * is guaranteed to produce results at least with the following relative accuracy
     */
    static double RELATIVE_ACCURACY = 1e-3;

    /**
     * The interpolation scheme used for calculating the log normalization factor of the emission probability
     * is guaranteed to produce results at least with the following absolute accuracy
     */
    static double ABSOLUTE_ACCURACY = 1e-3;

    /**
     * If read count is less than the following threshold, the Laplace approximation is not used.
     *
     * This choice of the threshold must coincide with the learning read count threshold since the learning
     * algorithm is based on the Laplace approximation.
     */
    private final int readCountThresholdPoissonSwitch;

    private final boolean applyPoissonToGaussianContinuityCorrection;

    /**
     * Various emission probability calculation models
     */
    public enum EmissionCalculationStrategy {
        /**
         * Fully Poisson; this mode is experimental since this model is incompatible with the EM
         * algorithm
         */
        POISSON,

        /**
         * Hybrid Poisson and Laplace+Poisson; this must be used for final posterior calling
         * (i.e. after learning). The Poisson approximation is used on low coverage targets
         * and Laplace approximation on Poisson is used on high coverage targets
         */
        HYBRID_POISSON_GAUSSIAN
    }

    /**
     * Public constructor.
     *
     * @param readCountThresholdPoissonSwitch read count threshold for switching between the exact and inexact
     *                                        Laplace-approximation-based calculation of the emission probability
     * @param applyPoissonToGaussianContinuityCorrection if true, continuity correction will be applied
     */
    public CoverageModelCopyRatioEmissionProbabilityCalculator(final int readCountThresholdPoissonSwitch,
                                                               final boolean applyPoissonToGaussianContinuityCorrection) {
        this.readCountThresholdPoissonSwitch = ParamUtils.isPositiveOrZero(readCountThresholdPoissonSwitch,
                "The read count threshold for switching to the exact Poisson emission calculation mode (i.e. no" +
                        " Laplace approximation) must be non-negative");
        this.applyPoissonToGaussianContinuityCorrection = applyPoissonToGaussianContinuityCorrection;
    }

    /**
     * Public constructor.
     *
     * @param readCountThresholdPoissonSwitch read count threshold for switching between the exact and inexact
     *                                        Laplace-approximation-based calculation of the emission probability
     */
    public CoverageModelCopyRatioEmissionProbabilityCalculator(final int readCountThresholdPoissonSwitch) {
        this.readCountThresholdPoissonSwitch = ParamUtils.isPositiveOrZero(readCountThresholdPoissonSwitch,
                "The read count threshold for switching to the exact Poisson emission calculation mode (i.e. no" +
                        " Laplace approximation) must be non-negative");
        this.applyPoissonToGaussianContinuityCorrection = true;
    }

    /**
     * Calculate the log emission probability. The parameter {@param target} neglected since
     * {@param emissionData} contains the necessary information.
     *
     * @param emissionData an instance of {@link CoverageModelCopyRatioEmissionData}
     * @param copyRatio copy ratio on which the emission probability is conditioned on
     * @param target target on which the emission probability is calculated (this parameter is not used)
     * @return emission probability
     */
    @Override
    public double logLikelihood(@Nonnull CoverageModelCopyRatioEmissionData emissionData, double copyRatio, @Nullable Target target) {
        Utils.nonNull(emissionData, "The emission data must be non-null");
        ParamUtils.isPositiveOrZero(copyRatio, "Copy ratio must be non-negative");

        final EmissionCalculationStrategy strategy =
                emissionData.getCopyRatioCallingMetadata().getEmissionCalculationStrategy();
        final boolean usePoisson = strategy.equals(EmissionCalculationStrategy.POISSON) ||
                emissionData.getReadCount() < readCountThresholdPoissonSwitch || copyRatio == 0;
        final double logLikelihood = usePoisson ? logLikelihoodPoisson(emissionData, copyRatio)
                                                : logLikelihoodLaplaceApproximation(emissionData, copyRatio);

        if (CHECK_FOR_NANS && Double.isNaN(logLikelihood)) {
            throw new RuntimeException("A NaN was produced while calculating the emission probability for" +
                    " emission data: " + emissionData.toString() + ", target: " + target +
                    ", copy ratio: " + copyRatio);
        }

        return logLikelihood;
    }

    /**
     * Calculate the log emission probability via Laplace approximation
     *
     * @param emissionData an instance of {@link CoverageModelCopyRatioEmissionData}
     * @param copyRatio copy ratio on which the emission probability is conditioned on
     * @return a double value
     */
    private double logLikelihoodLaplaceApproximation(@Nonnull final CoverageModelCopyRatioEmissionData emissionData, double copyRatio) {
        final double readDepth = emissionData.getCopyRatioCallingMetadata().getSampleCoverageDepth();
        final double err = emissionData.getMappingErrorProbability();
        final double copyRatioCorrection = err * FastMath.exp(-emissionData.getMu());
        final double mu = emissionData.getMu() + FastMath.log((1 - err) * readDepth * copyRatio +
                readDepth * copyRatioCorrection);
        final double psi = emissionData.getPsi();
        final double effectiveReadCount = applyPoissonToGaussianContinuityCorrection
                ? (double)emissionData.getReadCount() + 0.5
                : (double)emissionData.getReadCount();
        return getUnnormalizedLogProbabilityLaplaceApproximation(effectiveReadCount, mu, psi)
                - getLogProbabilityLaplaceApproximationNormalizationConstant(mu, psi);
    }

    /**
     * Calculate emission log probability directly using Poisson distribution
     *
     * TODO github/gatk-protected issue #854
     *
     * @implNote This is a naive implementation where the variance of log bias ($\Psi$) is
     * ignored altogether. This routine must be improved by performing a few-point numerical
     * integration of:
     *
     *      \int_{-\infty}^{+\infty} db Poisson(n | \lambda = d*c*exp(b) + eps*d)
     *                                                           * Normal(b | \mu, \Psi)
     *
     * @param emissionData an instance of {@link CoverageModelCopyRatioEmissionData}
     * @param copyRatio copy ratio on which the emission probability is conditioned on
     * @return a double value
     */
    private double logLikelihoodPoisson(@Nonnull CoverageModelCopyRatioEmissionData emissionData, double copyRatio) {
        final double multBias = FastMath.exp(emissionData.getMu());
        final double readDepth = emissionData.getCopyRatioCallingMetadata().getSampleCoverageDepth();
        final double err = emissionData.getMappingErrorProbability();
        final double poissonMean = readDepth * ((1 - err) * copyRatio * multBias + err);
        return new PoissonDistribution(null, poissonMean, PoissonDistribution.DEFAULT_EPSILON,
                PoissonDistribution.DEFAULT_MAX_ITERATIONS).logProbability(emissionData.getReadCount());
    }

    private static double getUnnormalizedLogProbabilityLaplaceApproximation(final double readCount, final double mu,
                                                                            final double psi) {
        final double variance = psi + 1.0 / readCount;
        final double logReadCount = FastMath.log(readCount);
        return - logReadCount - 0.5 * (LOG_2PI + FastMath.log(variance) + MathUtils.square(logReadCount - mu) / variance);
    }

    @VisibleForTesting
    static double getLogProbabilityLaplaceApproximationNormalizationConstant(final double mu, final double psi) {
        if (mu > MAXIMUM_TABULATED_MU) { /* the log norm factor is guaranteed to be smaller than 1e-4 */
            return 0.0;
        }
        final double muTrunc = FastMath.max(MINIMUM_TABULATED_MU, mu);
        final double psiTrunc = FastMath.max(MINIMUM_TABULATED_PSI, FastMath.min(MAXIMUM_TABULATED_PSI, psi));
        return logNormFactorSpline.value(muTrunc, psiTrunc);
    }

    /**
     * TODO github/gatk-protected issue #855 -- rewrite using org.broadinstitute.hellbender.utils.Utils.stream
     */
    private static double[] loadDoubleArrayTable(final InputStream inputStream) {
        final Scanner reader = new Scanner(inputStream);
        final List<Double> data = new ArrayList<>();
        while (reader.hasNextLine()) {
            data.add(Double.parseDouble(reader.nextLine()));
        }
        return data.stream().mapToDouble(d -> d).toArray();
    }

    /**
     * TODO github/gatk-protected issue #855 -- rewrite using org.broadinstitute.hellbender.utils.Utils.stream
     */
    private static double[][] loadDouble2DArrayTable(final InputStream inputStream) {
        final Scanner reader = new Scanner(inputStream);
        final List<double[]> data = new ArrayList<>();
        while (reader.hasNextLine()) {
            data.add(Arrays.stream(reader.nextLine().split("\t"))
                    .mapToDouble(Double::parseDouble).toArray());
        }
        final int rows = data.size();
        final int cols = data.get(0).length;
        final double[][] data2DArray = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                data2DArray[i][j] = data.get(i)[j];
            }
        }
        return data2DArray;
    }

    private static BicubicInterpolatingFunction getLogNormFactorInterpolatingFunction() {
        final Class<?> clazz = CoverageModelCopyRatioEmissionProbabilityCalculator.class;
        final double[] mu = loadDoubleArrayTable(clazz.getResourceAsStream(MU_TABLE_RESOURCE));
        final double[] psi = loadDoubleArrayTable(clazz.getResourceAsStream(PSI_TABLE_RESOURCE));
        final double[][] logNorm = loadDouble2DArrayTable(clazz.getResourceAsStream(LOG_NORM_TABLE_RESOURCE));
        final BicubicInterpolator interp = new BicubicInterpolator();
        return interp.interpolate(mu, psi, logNorm);
    }
}
