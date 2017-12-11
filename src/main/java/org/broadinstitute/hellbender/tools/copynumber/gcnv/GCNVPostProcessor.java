package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nullable;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

/**
 * Helper class for single sample gCNV postprocessing
 */
public class GCNVPostProcessor {

    /**
     * VCF header keys
     */
    public static final String CN_MAP = "CNMAP";
    public static final String CN_MEAN = "CNMEAN";
    public static final String CN_STD = "CNSTD";
    public static final String CNQ = "CNQ";

    private final List<Allele> alleles;
    private final IntegerCopyNumberStateCollection integerCopyNumberStateCollection;
    private final String sampleName;
    private final VariantContextWriter outputWriter;

    GCNVPostProcessor(final VariantContextWriter outputWriter,
                      final IntegerCopyNumberStateCollection integerCopyNumberStateCollection,
                      final String sampleName) {
        this.outputWriter = Utils.nonNull(outputWriter);
        this.integerCopyNumberStateCollection = Utils.nonNull(integerCopyNumberStateCollection);
        this.alleles = integerCopyNumberStateCollection.getAlleles();
        this.sampleName = sampleName;
    }

    /**
     *
     */
    public void composeVariantContextAndWrite(@Nullable final String commandLine) {
        //TODO pass a list of intervals and add progress meter
        //ProgressMeter progressMeter = new ProgressMeter(1.0);
        //progressMeter.start();
        outputWriter.writeHeader(composeHeader(commandLine));
    }

    public void writeChunkedVariantContext(final List<CopyNumberPosteriorLocatableRecord> copyNumberPosteriorLocatableRecordsChunk,
                                           final String variantPrefix) {
        for (CopyNumberPosteriorLocatableRecord posteriorRecord: copyNumberPosteriorLocatableRecordsChunk) {
            final VariantContext variantContext = composeVariantContext(posteriorRecord, variantPrefix);
            outputWriter.add(variantContext);
        }
    }

    /**
     * Composes the VCF header
     */
    private VCFHeader composeHeader(@Nullable final String commandLine) {
        final VCFHeader result = new VCFHeader(Collections.emptySet(), Arrays.asList(sampleName));

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        /* add command line */
        if (commandLine != null) {
            result.addMetaDataLine(new VCFHeaderLine("command", commandLine));
        }

        /* header lines related to formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(CN_MAP, 1,
                VCFHeaderLineType.Integer, "Copy number maximum a posteriori"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CN_MEAN, 1,
                VCFHeaderLineType.Float, "Copy number posterior mean"));
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                VCFHeaderLineType.Integer, "End coordinate of this variant"));

        return result;
    }

    /**
     *
     * @param copyNumberPosteriorLocatableRecord a posterior record to genotype
     * @param variantPrefix a variant prefix
     * @return composed variant context
     */
    protected VariantContext composeVariantContext(final CopyNumberPosteriorLocatableRecord copyNumberPosteriorLocatableRecord,
                                                 final String variantPrefix) {
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.alleles(alleles);
        variantContextBuilder.chr(copyNumberPosteriorLocatableRecord.getContig());
        variantContextBuilder.start(copyNumberPosteriorLocatableRecord.getStart());
        variantContextBuilder.stop(copyNumberPosteriorLocatableRecord.getEnd());
        variantContextBuilder.id(String.format(variantPrefix + "_%s_%d_%d",
                copyNumberPosteriorLocatableRecord.getContig(),
                copyNumberPosteriorLocatableRecord.getStart(),
                copyNumberPosteriorLocatableRecord.getEnd()));
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final int copyNumberMAP = calculateMAPCopyNumber(copyNumberPosteriorLocatableRecord);
        genotypeBuilder.attribute(CN_MAP, copyNumberMAP);
        final Genotype genotype = genotypeBuilder.make();
        //TODO add variant context fields to genotype

        //Add allele information to the variant context
        variantContextBuilder.attribute(VCFConstants.END_KEY, copyNumberPosteriorLocatableRecord.getEnd());

        variantContextBuilder.genotypes(genotype);
        return variantContextBuilder.make();
    }

    private int calculateMAPCopyNumber(CopyNumberPosteriorLocatableRecord copyNumberPosteriorLocatableRecord) {
        Optional<IntegerCopyNumberState> copyNumberStateMAP = integerCopyNumberStateCollection.getCopyNumberStates()
                .stream().max((a1, a2) -> Double.compare(copyNumberPosteriorLocatableRecord.getCopyNumberPosteriors().getCopyNumberPosterior(a1),
                copyNumberPosteriorLocatableRecord.getCopyNumberPosteriors().getCopyNumberPosterior(a2)));
        return copyNumberStateMAP.get().getCopyNumber();
    }


}
