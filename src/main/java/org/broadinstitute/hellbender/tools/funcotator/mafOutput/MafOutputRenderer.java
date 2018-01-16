package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.Funcotator;
import org.broadinstitute.hellbender.tools.funcotator.OutputRenderer;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * A Funcotator output renderer for writing to MAF files.
 * Complies with version 2.4 of the MAF standard found here:
 * https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification+-+v2.4
 * Created by jonn on 12/5/17.
 */
public class MafOutputRenderer extends OutputRenderer {

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    /**
     * Value to insert into columns where no annotation exists.
     */
    protected static final String UNKNOWN_VALUE_STRING = "__UNKNOWN__";

    /**
     * Default set of columns to include in this {@link MafOutputRenderer}.
     * Order of the columns is preserved by the {@link LinkedHashMap}, while still being able to access each field via
     * the associated key.
     */
    private static final LinkedHashMap<String, String> defaultMap = new LinkedHashMap<>();

    /**
     * Delimiter for fields in the output MAF file.
     */
    private static final String FIELD_DELIMITER = "\t";

    //==================================================================================================================
    // Private Members:

    /**
     * {@link List} of the {@link DataSourceFuncotationFactory} objects that are being used in this run of {@link Funcotator}.
     */
    private final List<DataSourceFuncotationFactory> dataSourceFactories;

    /**
     * {@link Path} to the resulting MAF output file for this {@link MafOutputRenderer}.
     */
    private final Path outputFilePath;

    //==================================================================================================================
    // Constructors:

    public MafOutputRenderer(final Path outputFilePath,
                             final List<DataSourceFuncotationFactory> dataSources,
                             final LinkedHashMap<String, String> unaccountedForDefaultAnnotations,
                             final LinkedHashMap<String, String> unaccountedForOverrideAnnotations) {

        this.outputFilePath = outputFilePath;

        dataSourceFactories = dataSources;

        // Merge the annotations into our manualAnnotations:
        manualAnnotations = new LinkedHashMap<>();
        manualAnnotations.putAll(unaccountedForDefaultAnnotations);
        manualAnnotations.putAll(unaccountedForOverrideAnnotations);

        // Cache the manual annotation string so we can pass it easily into any Funcotations:
        manualAnnotationSerializedString = (manualAnnotations.size() != 0 ? String.join( FIELD_DELIMITER, manualAnnotations.values() ) + FIELD_DELIMITER : "");

        initializeDefaultMapWithKeys();
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public void open() {

    }

    @Override
    public void close() {

    }

    @Override
    public void write(final VariantContext variant, final List<Funcotation> funcotations) {

        //TODO: For this to work correctly you'll need to group funcotations by reference / alternate allele pairs.

        // First separate the funcotations into Gencode and Other:
        final List<GencodeFuncotation> gencodeFuncotations = new ArrayList<>(funcotations.size() / 2);
        final List<Funcotation> otherFuncotations = new ArrayList<>(funcotations.size() / 2);
        for ( final Funcotation funcotation : funcotations ) {
            if ( GencodeFuncotation.class.isAssignableFrom(funcotation.getClass()) ) {
                gencodeFuncotations.add((GencodeFuncotation) funcotation);
            }
            else {
                otherFuncotations.add(funcotation);
            }
        }

        // Create our output map:
        final LinkedHashMap<String, Object> outputMap = new LinkedHashMap<>(defaultMap);

        // Go through Gencode funcotations as they're the basis for each line:
        for ( final GencodeFuncotation gencodeFuncotation : gencodeFuncotations ) {
            outputMap.put("Hugo_symbol", gencodeFuncotation.getHugoSymbol());
            outputMap.put("NCBI_Build", gencodeFuncotation.getNcbiBuild());
            outputMap.put("Chromosome", gencodeFuncotation.getChromosome());
            outputMap.put("Start_Position", gencodeFuncotation.getStart());
            outputMap.put("End_Position", gencodeFuncotation.getEnd());
            outputMap.put("Strand", gencodeFuncotation.getTranscriptStrand());
            outputMap.put("VariantClassification", gencodeFuncotation.getVariantClassification());
            outputMap.put("Variant_Type", gencodeFuncotation.getVariantType());
            outputMap.put("Reference_Allele", gencodeFuncotation.getRefAllele());
            outputMap.put("Tumor_Seq_Allele1", gencodeFuncotation.getTumorSeqAllele1());
            outputMap.put("Tumor_Seq_Allele2", gencodeFuncotation.getTumorSeqAllele2());
        }

        for ( final Funcotation funcotation : otherFuncotations ) {
            //TODO: For this to work you'll need to make annotations work as maps ala issue #3919
        }

    }


    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    /**
     * Initializes {@link MafOutputRenderer#defaultMap} with the default keys for the columns in a MAF file.
     */
    void initializeDefaultMapWithKeys() {
        defaultMap.put("Hugo_Symbol",                   UNKNOWN_VALUE_STRING );
        defaultMap.put("Entrez_Gene_Id",                UNKNOWN_VALUE_STRING );
        defaultMap.put("Center",                        UNKNOWN_VALUE_STRING );
        defaultMap.put("NCBI_Build",                    UNKNOWN_VALUE_STRING );
        defaultMap.put("Chromosome",                    UNKNOWN_VALUE_STRING );
        defaultMap.put("Start_Position",                UNKNOWN_VALUE_STRING );
        defaultMap.put("End_Position",                  UNKNOWN_VALUE_STRING );
        defaultMap.put("Strand",                        UNKNOWN_VALUE_STRING );
        defaultMap.put("Variant_Classification",        UNKNOWN_VALUE_STRING );
        defaultMap.put("Variant_Type",                  UNKNOWN_VALUE_STRING );
        defaultMap.put("Reference_Allele",              UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Seq_Allele1",             UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Seq_Allele2",             UNKNOWN_VALUE_STRING );
        defaultMap.put("dbSNP_RS",                      UNKNOWN_VALUE_STRING );
        defaultMap.put("dbSNP_Val_Status",              UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Sample_Barcode",          UNKNOWN_VALUE_STRING );
        defaultMap.put("Matched_Norm_Sample_Barcode",   UNKNOWN_VALUE_STRING );
        defaultMap.put("Match_Norm_Seq_Allele1",        UNKNOWN_VALUE_STRING );
        defaultMap.put("Match_Norm_Seq_Allele2",        UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Validation_Allele1",      UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Validation_Allele2",      UNKNOWN_VALUE_STRING );
        defaultMap.put("Match_Norm_Validation_Allele1", UNKNOWN_VALUE_STRING );
        defaultMap.put("Match_Norm_Validation_Allele2", UNKNOWN_VALUE_STRING );
        defaultMap.put("Verification_Status",           UNKNOWN_VALUE_STRING );
        defaultMap.put("Validation_Status",             UNKNOWN_VALUE_STRING );
        defaultMap.put("Mutation_Status",               UNKNOWN_VALUE_STRING );
        defaultMap.put("Sequencing_Phase",              UNKNOWN_VALUE_STRING );
        defaultMap.put("Sequence_Source",               UNKNOWN_VALUE_STRING );
        defaultMap.put("Validation_Method",             UNKNOWN_VALUE_STRING );
        defaultMap.put("Score",                         UNKNOWN_VALUE_STRING );
        defaultMap.put("BAM_File",                      UNKNOWN_VALUE_STRING );
        defaultMap.put("Sequencer",                     UNKNOWN_VALUE_STRING );
        defaultMap.put("Tumor_Sample_UUID",             UNKNOWN_VALUE_STRING );
        defaultMap.put("Matched_Norm_Sample_UUID",      UNKNOWN_VALUE_STRING );
    }

    //==================================================================================================================
    // Helper Data Types:

    //------------------------------------------------------------------------------------------------------------------
    // Required Columns:

// 1	Hugo_Symbol
// 2	Entrez_Gene_Id
// 3	Center
// 4	NCBI_Build
// 5	Chromosome
// 6	Start_Position
// 7	End_Position
// 8	Strand
// 9	Variant_Classification
//10	Variant_Type
//11	Reference_Allele
//12	Tumor_Seq_Allele1
//13	Tumor_Seq_Allele2
//14	dbSNP_RS
//15	dbSNP_Val_Status
//16	Tumor_Sample_Barcode
//17	Matched_Norm_Sample_Barcode
//18	Match_Norm_Seq_Allele1
//19	Match_Norm_Seq_Allele2
//20	Tumor_Validation_Allele1
//21	Tumor_Validation_Allele2
//22	Match_Norm_Validation_Allele1
//23	Match_Norm_Validation_Allele2
//24	Verification_Status
//25	Validation_Status
//26	Mutation_Status
//27	Sequencing_Phase
//28	Sequence_Source
//29	Validation_Method
//30	Score
//31	BAM_File
//32	Sequencer
//33    Tumor_Sample_UUID
//34    Matched_Norm_Sample_UUID

}
