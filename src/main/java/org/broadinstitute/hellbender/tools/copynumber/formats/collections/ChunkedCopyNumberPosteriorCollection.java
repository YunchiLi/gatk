package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.gcnv.CopyNumberPosteriorRecord;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberStateCollection;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Collection of copy number posteriors for an individual chunk containing a subset of intervals considered in the analysis
 */
public class ChunkedCopyNumberPosteriorCollection extends AbstractSampleRecordCollection<CopyNumberPosteriorRecord> {

    public ChunkedCopyNumberPosteriorCollection(final File inputFile,
                                                final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        super(inputFile,
                Utils.nonNull(integerCopyNumberStateCollection).getTableColumnCollection(),
                getPosteriorRecordFromDataLineDecoder(integerCopyNumberStateCollection),
                getPosteriorRecordToDataLineEncoder(integerCopyNumberStateCollection));
    }

    private static Function<DataLine, CopyNumberPosteriorRecord> getPosteriorRecordFromDataLineDecoder(final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        return dataLine -> {
            final Map<IntegerCopyNumberState, Double> copyNumberStateDoubleMap = new HashMap<>();
            for(int i = 0; i < integerCopyNumberStateCollection.size(); i++) {
                copyNumberStateDoubleMap.putIfAbsent(integerCopyNumberStateCollection.get(i), dataLine.getDouble(i));
            }
            final CopyNumberPosteriorRecord record = new CopyNumberPosteriorRecord(copyNumberStateDoubleMap);
            return record;
        };
    }

    private static BiConsumer<CopyNumberPosteriorRecord, DataLine> getPosteriorRecordToDataLineEncoder(final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        //TODO
        return (copyNumberPosteriorRecord, dataLine) -> {
            return;
        };
    }

}
