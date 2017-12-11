package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;

/**
 * A copy number posterior record
 */
public class CopyNumberPosteriorRecord {
    private final Map<IntegerCopyNumberState, Double> copyNumberStatePosteriors;

    public CopyNumberPosteriorRecord(final Map<IntegerCopyNumberState, Double> copyNumberStatePosteriors) {
        this.copyNumberStatePosteriors = Utils.nonNull(copyNumberStatePosteriors);
    }

    /**
     * Get the posterior probability for a given copy number state
     */
    public Double getCopyNumberPosterior(IntegerCopyNumberState integerCopyNumberState) {
        return copyNumberStatePosteriors.get(integerCopyNumberState);
    }

}
