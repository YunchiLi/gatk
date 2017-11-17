package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.avro.generic.GenericData;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.VariantAnnotationArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.utils.ClassUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.reflections.ReflectionUtils;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * The class responsible for computing annotations for variants.
 * Annotations are auto-discovered - ie, any class that extends {@link VariantAnnotation} and
 * lives in this package is treated as an annotation and the engine will attempt to create instances of it
 * by calling the non-arg constructor (loading will fail if there is no no-arg constructor).
 */
public final class VariantAnnotatorEngine {
    private final List<InfoFieldAnnotation> infoAnnotations;
    private final List<GenotypeAnnotation> genotypeAnnotations;
    private Set<String> reducibleKeys;
    private List<VAExpression> expressions = new ArrayList<>();

    private final VariantOverlapAnnotator variantOverlapAnnotator;
    private Boolean expressionAlleleConcordance;
    private final Boolean useRawAnnotations;

//    private VariantAnnotatorEngine(final AnnotationManager annots,
//                                   final FeatureInput<VariantContext> dbSNPInput,
//                                   final List<FeatureInput<VariantContext>> featureInputs,
//                                   final boolean useRaw){
//        infoAnnotations = annots.createInfoFieldAnnotations();
//        genotypeAnnotations = annots.createGenotypeAnnotations();
//        variantOverlapAnnotator = initializeOverlapAnnotator(dbSNPInput, featureInputs);
//        reducibleKeys = new HashSet<>();
//        useRawAnnotations = useRaw;
//        for (InfoFieldAnnotation annot : infoAnnotations) {
//            if (annot instanceof ReducibleAnnotation) {
//                reducibleKeys.add(((ReducibleAnnotation) annot).getRawKeyName());
//            }
//        }
//    }

    /**
     * TODO comment this out
     */
    public VariantAnnotatorEngine(final List<Annotation> annotationList,
                                   final FeatureInput<VariantContext> dbSNPInput,
                                   final List<FeatureInput<VariantContext>> featureInputs,
                                   final boolean useRaw){
        infoAnnotations = new ArrayList<>();
        genotypeAnnotations = new ArrayList<>();
        for (Annotation annot : annotationList) {
            if (annot instanceof InfoFieldAnnotation) {
                infoAnnotations.add((InfoFieldAnnotation) annot);
            }
            if (annot instanceof GenotypeAnnotation) {
                genotypeAnnotations.add((GenotypeAnnotation) annot);
            }
        }
        variantOverlapAnnotator = initializeOverlapAnnotator(dbSNPInput, featureInputs);
        reducibleKeys = new HashSet<>();
        useRawAnnotations = useRaw;
        for (InfoFieldAnnotation annot : infoAnnotations) {
            if (annot instanceof ReducibleAnnotation) {
                reducibleKeys.add(((ReducibleAnnotation) annot).getRawKeyName());
            }
        }
    }


    /**
     * Makes the engine for all known annotation types (minus the excluded ones).
     * @param annotationsToExclude list of annotations to exclude (pass an empty list to indicate that there are no exclusions)
     * @param dbSNPInput input for variants from a known set from DbSNP or null if not provided.
     *                   The annotation engine will mark variants overlapping anything in this set using {@link htsjdk.variant.vcf.VCFConstants#DBSNP_KEY}.
     * @param comparisonFeatureInputs list of inputs with known variants.
     *                   The annotation engine will mark variants overlapping anything in those sets using the name given by {@link FeatureInput#getName()}.
     *                   Note: the DBSNP FeatureInput should be passed in separately, and not as part of this List - an GATKException will be thrown otherwise.
     *                   Note: there are no non-DBSNP comparison FeatureInputs an empty List should be passed in here, rather than null.
     */
    public static VariantAnnotatorEngine ofAllMinusExcluded(final List<String> annotationsToExclude,
                                                            final FeatureInput<VariantContext> dbSNPInput,
                                                            final List<FeatureInput<VariantContext>> comparisonFeatureInputs,
                                                            final Boolean useRawAnnotations) {
        Utils.nonNull(annotationsToExclude, "annotationsToExclude is null");
        Utils.nonNull(comparisonFeatureInputs, "comparisonFeatureInputs is null");
        return null;//new VariantAnnotatorEngine(AnnotationManager.ofAllMinusExcluded(annotationsToExclude), dbSNPInput, comparisonFeatureInputs, useRawAnnotations);
    }

    /**
     * Makes the engine for given annotation types and groups (minus the excluded ones).
     * @param annotationGroupsToUse list of annotations groups to include
     * @param annotationsToUse     list of of annotations to include
     * @param annotationsToExclude list of annotations to exclude
     * @param dbSNPInput input for variants from a known set from DbSNP or null if not provided.
     *                   The annotation engine will mark variants overlapping anything in this set using {@link htsjdk.variant.vcf.VCFConstants#DBSNP_KEY}.
     * @param comparisonFeatureInputs list of inputs with known variants.
     *                   The annotation engine will mark variants overlapping anything in those sets using the name given by {@link FeatureInput#getName()}.
     *                   Note: the DBSNP FeatureInput should be passed in separately, and not as part of this List - an GATKException will be thrown otherwise.
     *                   Note: there are no non-DBSNP comparison FeatureInputs an empty List should be passed in here, rather than null.
     */
    public static VariantAnnotatorEngine ofSelectedMinusExcluded(final List<String> annotationGroupsToUse,
                                                                 final List<String> annotationsToUse,
                                                                 final List<String> annotationsToExclude,
                                                                 final FeatureInput<VariantContext> dbSNPInput,
                                                                 final List<FeatureInput<VariantContext>> comparisonFeatureInputs) {
        Utils.nonNull(annotationGroupsToUse, "annotationGroupsToUse is null");
        Utils.nonNull(annotationsToUse, "annotationsToUse is null");
        Utils.nonNull(annotationsToExclude, "annotationsToExclude is null");
        Utils.nonNull(comparisonFeatureInputs, "comparisonFeatureInputs is null");
        return null;//new VariantAnnotatorEngine(AnnotationManager.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude), dbSNPInput, comparisonFeatureInputs, false);
    }

    /**
     * Makes the engine for given annotation types and groups (minus the excluded ones).
     * @param annotationGroupsToUse list of annotations groups to include
     * @param annotationsToUse     list of of annotations to include
     * @param annotationsToExclude list of annotations to exclude
     * @param dbSNPInput input for variants from a known set from DbSNP or null if not provided.
     *                   The annotation engine will mark variants overlapping anything in this set using {@link htsjdk.variant.vcf.VCFConstants#DBSNP_KEY}.
     * @param comparisonFeatureInputs list of inputs with known variants.
     *                   The annotation engine will mark variants overlapping anything in those sets using the name given by {@link FeatureInput#getName()}.
     *                   Note: the DBSNP FeatureInput should be passed in separately, and not as part of this List - an GATKException will be thrown otherwise.
     *                   Note: there are no non-DBSNP comparison FeatureInputs an empty List should be passed in here, rather than null.
     * @param useRawAnnotations  Specify whether the annotation engine will attempt to output raw annotations for reducible annotations
     */
    public static VariantAnnotatorEngine ofSelectedMinusExcluded(final List<String> annotationGroupsToUse,
                                                                 final List<String> annotationsToUse,
                                                                 final List<String> annotationsToExclude,
                                                                 final FeatureInput<VariantContext> dbSNPInput,
                                                                 final List<FeatureInput<VariantContext>> comparisonFeatureInputs,
                                                                 final Boolean useRawAnnotations) {
        Utils.nonNull(annotationGroupsToUse, "annotationGroupsToUse is null");
        Utils.nonNull(annotationsToUse, "annotationsToUse is null");
        Utils.nonNull(annotationsToExclude, "annotationsToExclude is null");
        Utils.nonNull(comparisonFeatureInputs, "comparisonFeatureInputs is null");
        return null;//VariantAnnotatorEngine(AnnotationManager.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude), dbSNPInput, comparisonFeatureInputs, useRawAnnotations);
    }

    /**
     * An overload of {@link org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine#ofSelectedMinusExcluded ofSelectedMinusExcluded}
     * except that it accepts a {@link VariantAnnotationArgumentCollection} as input.
     * @param argumentCollection            VariantAnnotationArgumentCollection containing requested annotations.
     * @param dbSNPInput                    input for variants from a known set from DbSNP or null if not provided.
     *                   The annotation engine will mark variants overlapping anything in this set using {@link htsjdk.variant.vcf.VCFConstants#DBSNP_KEY}.
     * @param comparisonFeatureInputs list of inputs with known variants.
     *                   The annotation engine will mark variants overlapping anything in those sets using the name given by {@link FeatureInput#getName()}.
     *                   Note: the DBSNP FeatureInput should be passed in separately, and not as part of this List - an GATKException will be thrown otherwise.
     *                   Note: there are no non-DBSNP comparison FeatureInputs an empty List should be passed in here, rather than null.
     * @return a VariantAnnotatorEngine initialized with the requested annotations
     */
    public static VariantAnnotatorEngine ofSelectedMinusExcluded(final VariantAnnotationArgumentCollection argumentCollection,
                                                                 final FeatureInput<VariantContext> dbSNPInput,
                                                                 final List<FeatureInput<VariantContext>> comparisonFeatureInputs,
                                                                 final Boolean useRawAnnotations) {
        return ofSelectedMinusExcluded(argumentCollection.annotationGroupsToUse,
                argumentCollection.annotationsToUse,
                argumentCollection.annotationsToExclude,
                dbSNPInput, comparisonFeatureInputs, useRawAnnotations);
    }
    private VariantOverlapAnnotator initializeOverlapAnnotator(final FeatureInput<VariantContext> dbSNPInput, final List<FeatureInput<VariantContext>> featureInputs) {
        final Map<FeatureInput<VariantContext>, String> overlaps = new LinkedHashMap<>();
        for ( final FeatureInput<VariantContext> fi : featureInputs) {
            overlaps.put(fi, fi.getName());
        }
        if (overlaps.values().contains(VCFConstants.DBSNP_KEY)){
            throw new GATKException("The map of overlaps must not contain " + VCFConstants.DBSNP_KEY);
        }
        if (dbSNPInput != null) {
            overlaps.put(dbSNPInput, VCFConstants.DBSNP_KEY); // add overlap detection with DBSNP by default
        }

        return new VariantOverlapAnnotator(dbSNPInput, overlaps);
    }

    //TODO this is a hack and should not be tried at home kids
    public void removeInfoAnnotation(Class<? extends InfoFieldAnnotation> toRemove) {
        for (InfoFieldAnnotation annotation : infoAnnotations) {
            if (annotation.getClass() == toRemove) {
                infoAnnotations.remove(annotation);
                return;
            }
        }
    }

    /**
     * Returns the list of genotype annotations that will be applied.
     * Note: The returned list is unmodifiable.
     */
    public List<GenotypeAnnotation> getGenotypeAnnotations() {
        return Collections.unmodifiableList(genotypeAnnotations);
    }

    /**
     * Returns the list of info annotations that will be applied.
     * Note: The returned list is unmodifiable.
     */
    public List<InfoFieldAnnotation> getInfoAnnotations() {
        return Collections.unmodifiableList(infoAnnotations);
    }

    /**
     * Returns the set of descriptions to be added to the VCFHeader line (for all annotations in this engine).
     */
    public Set<VCFHeaderLine> getVCFAnnotationDescriptions() {
        return getVCFAnnotationDescriptions(false);
    }

    /**
     * Returns the set of descriptions to be added to the VCFHeader line (for all annotations in this engine).
     * @param useRaw Whether to prefer reducible annotation raw key descriptions over their normal descriptions
     */
    public Set<VCFHeaderLine> getVCFAnnotationDescriptions(boolean useRaw) {
        final Set<VCFHeaderLine> descriptions = new LinkedHashSet<>();

        for ( final InfoFieldAnnotation annotation : infoAnnotations) {
            if (annotation instanceof ReducibleAnnotation && useRaw) {
                descriptions.addAll(((ReducibleAnnotation)annotation).getRawDescriptions());
            } else {
                descriptions.addAll(annotation.getDescriptions());
            }
        }
        for ( final GenotypeAnnotation annotation : genotypeAnnotations) {
            descriptions.addAll(annotation.getDescriptions());
        }
        for ( final String db : variantOverlapAnnotator.getOverlapNames() ) {
            if ( VCFStandardHeaderLines.getInfoLine(db, false) != null ) {
                descriptions.add(VCFStandardHeaderLines.getInfoLine(db));
            } else {
                descriptions.add(new VCFInfoHeaderLine(db, 0, VCFHeaderLineType.Flag, db + " Membership"));
            }
        }

        Utils.validate(!descriptions.contains(null), "getVCFAnnotationDescriptions should not contain null. This error is likely due to an incorrect implementation of getDescriptions() in one or more of the annotation classes");
        return descriptions;
    }


    /**
     * Combine (raw) data for reducible annotations (those that use raw data in gVCFs)
     * Mutates annotationMap by removing the annotations that were combined
     *
     * Additionally, will combine other annotations by parsing them as numbers and reducing them
     * down to the
     * @param allelesList   the list of merged alleles across all variants being combined
     * @param annotationMap attributes of merged variant contexts -- is modifying by removing successfully combined annotations
     * @return  a map containing the keys and raw values for the combined annotations
     */
    @SuppressWarnings({"unchecked"})
    public Map<String, Object> combineAnnotations(final List<Allele> allelesList, Map<String, List<?>> annotationMap) {
        Map<String, Object> combinedAnnotations = new HashMap<>();

        // go through all the requested reducible info annotationTypes
        for (final InfoFieldAnnotation annotationType : infoAnnotations) {
            if (annotationType instanceof ReducibleAnnotation) {
                ReducibleAnnotation currentASannotation = (ReducibleAnnotation) annotationType;
                if (annotationMap.containsKey(currentASannotation.getRawKeyName())) {
                    final List<ReducibleAnnotationData<?>> annotationValue = (List<ReducibleAnnotationData<?>>) annotationMap.get(currentASannotation.getRawKeyName());
                    final Map<String, Object> annotationsFromCurrentType = currentASannotation.combineRawData(allelesList, annotationValue);
                    combinedAnnotations.putAll(annotationsFromCurrentType);
                    //remove the combined annotations so that the next method only processes the non-reducible ones
                    annotationMap.remove(currentASannotation.getRawKeyName());
                }
            }
        }
        return combinedAnnotations;
    }

    /**
     * Finalize reducible annotations (those that use raw data in gVCFs)
     * @param vc    the merged VC with the final set of alleles, possibly subset to the number of maxAltAlleles for genotyping
     * @param originalVC    the merged but non-subset VC that contains the full list of merged alleles
     * @return  a VariantContext with the final annotation values for reducible annotations
     */
    public VariantContext finalizeAnnotations(VariantContext vc, VariantContext originalVC) {
        final Map<String, Object> variantAnnotations = new LinkedHashMap<>(vc.getAttributes());

        // go through all the requested info annotationTypes
        for (final InfoFieldAnnotation annotationType : infoAnnotations) {
            if (annotationType instanceof ReducibleAnnotation) {

                ReducibleAnnotation currentASannotation = (ReducibleAnnotation) annotationType;

                final Map<String, Object> annotationsFromCurrentType = currentASannotation.finalizeRawData(vc, originalVC);
                if (annotationsFromCurrentType != null) {
                    variantAnnotations.putAll(annotationsFromCurrentType);
                }
                //clean up raw annotation data after annotations are finalized
                variantAnnotations.remove(currentASannotation.getRawKeyName());
            }
        }

        // generate a new annotated VC
        final VariantContextBuilder builder = new VariantContextBuilder(vc).attributes(variantAnnotations);

        // annotate genotypes, creating another new VC in the process
        final VariantContext annotated = builder.make();
        return annotated;
    }

    /**
     * Annotates the given variant context - adds all annotations that satisfy the predicate.
     * @param vc the variant context to annotate
     * @param features context containing the features that overlap the given variant
     * @param ref the reference context of the variant to annotate or null if there is none
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample. May be null
     * @param addAnnot function that indicates if the given annotation type should be added to the variant
     *
     */
    public VariantContext annotateContext(final VariantContext vc,
                                          final FeatureContext features,
                                          final ReferenceContext ref,
                                          final ReadLikelihoods<Allele> likelihoods,
                                          final Predicate<VariantAnnotation> addAnnot) {
        Utils.nonNull(vc, "vc cannot be null");
        Utils.nonNull(features, "features cannot be null");
        Utils.nonNull(addAnnot, "addAnnot cannot be null");

        // annotate genotypes, creating another new VC in the process
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        builder.genotypes(annotateGenotypes(ref, vc, likelihoods, addAnnot));
        final VariantContext newGenotypeAnnotatedVC = builder.make();

        final Map<String, Object> infoAnnotMap = new LinkedHashMap<>(newGenotypeAnnotatedVC.getAttributes());
        annotateExpressions(vc, features, ref, infoAnnotMap);

        for ( final InfoFieldAnnotation annotationType : this.infoAnnotations) {
            if (addAnnot.test(annotationType)){
                final Map<String, Object> annotationsFromCurrentType;
                if (useRawAnnotations && annotationType instanceof ReducibleAnnotation) {
                    annotationsFromCurrentType = ((ReducibleAnnotation) annotationType).annotateRawData(ref, newGenotypeAnnotatedVC, likelihoods);
                } else {
                    annotationsFromCurrentType = annotationType.annotate(ref, newGenotypeAnnotatedVC, likelihoods);
                }
                if ( annotationsFromCurrentType != null ) {
                    infoAnnotMap.putAll(annotationsFromCurrentType);
                }
            }
        }

        // create a new VC with info and genotype annotations
        final VariantContext annotated = builder.attributes(infoAnnotMap).make();

        // annotate db occurrences
        return variantOverlapAnnotator.annotateOverlaps(features, variantOverlapAnnotator.annotateRsID(features, annotated));
    }

    private GenotypesContext annotateGenotypes(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final ReadLikelihoods<Allele> likelihoods,
                                               final Predicate<VariantAnnotation> addAnnot) {
        if ( genotypeAnnotations.isEmpty() ) {
            return vc.getGenotypes();
        }

        final GenotypesContext genotypes = GenotypesContext.create(vc.getNSamples());
        for ( final Genotype genotype : vc.getGenotypes() ) {
            final GenotypeBuilder gb = new GenotypeBuilder(genotype);
            for ( final GenotypeAnnotation annotation : genotypeAnnotations) {
                if (addAnnot.test(annotation)) {
                    annotation.annotate(ref, vc, genotype, gb, likelihoods);
                }
            }
            genotypes.add(gb.make());
        }

        return genotypes;
    }

    /**
     * Method which checks if a key is a raw key of the requested reducible annotations
     * @param key annotation key to check
     * @return true if the key is the raw key for a requested annotation
     */
    public boolean isRequestedReducibleRawKey(String key) {
        return reducibleKeys.contains(key);
    }

    private static final class AnnotationManager {

        private final List<String> annotationGroupsToUse;
        private final List<String> annotationsToUse;
        private final List<String> annotationsToExclude;

        private AnnotationManager(final List<String> annotationGroupsToUse, final List<String> annotationsToUse, final List<String> annotationsToExclude){
            this.annotationGroupsToUse = annotationGroupsToUse;
            this.annotationsToUse = annotationsToUse;
            this.annotationsToExclude = annotationsToExclude;

            final Set<String> allAnnotationNames = new LinkedHashSet<>(AnnotationManager.getAllAnnotationNames());
            final Set<String> unknownAnnots = Sets.difference(new LinkedHashSet<>(annotationsToUse), allAnnotationNames);
            assertAnnotationExists(unknownAnnots);

            final Set<String> unknownAnnotsExclude = Sets.difference(new LinkedHashSet<>(annotationsToExclude), allAnnotationNames);
            assertAnnotationExists(unknownAnnotsExclude);

            final Set<String> unknownGroups =  Sets.difference(new LinkedHashSet<>(annotationGroupsToUse), new LinkedHashSet<>(AnnotationManager.getAllAnnotationGroupNames()));
            if (!unknownGroups.isEmpty()){
                throw new CommandLineException.BadArgumentValue("group", "Unknown annotation group in " + unknownGroups + ". Known groups are " + AnnotationManager.getAllAnnotationGroupNames());
            }
        }

        private void assertAnnotationExists(final Set<String> missingAnnots){
            if (!missingAnnots.isEmpty()){
                throw new CommandLineException.BadArgumentValue("annotation", "Annotation " + missingAnnots + " not found; please check that you have specified the name correctly");
            }
        }

        /**
         * An annotation will be included only when:
         * - it is in one of the annotation groups or
         * - it is listed explicitly
         * - and it is not excluded explicitly.
         */
        static AnnotationManager ofSelectedMinusExcluded(final List<String> annotationGroupsToUse, final List<String> annotationsToUse, final List<String> annotationsToExclude){
            final List<String> groups = new ArrayList<>(annotationGroupsToUse);//make copy
            final List<String> annots = new ArrayList<>(annotationsToUse);//make copy
            final List<String> excludes = new ArrayList<>(annotationsToExclude);//make copy
            return new AnnotationManager(groups, annots, excludes);
        }

        /**
         * An annotation will be included only when it is not excluded explicitly.
         */
        static AnnotationManager ofAllMinusExcluded(final List<String> annotationsToExclude){
            final List<String> groups = getAllAnnotationGroupNames();
            final List<String> annots = getAllAnnotationNames();
            return new AnnotationManager(groups, annots, annotationsToExclude);
        }

        private static List<String> getAllAnnotationNames() {
            final Set<VariantAnnotation> union = Sets.union(new LinkedHashSet<>(makeAllGenotypeAnnotations()), new LinkedHashSet<>(AnnotationManager.makeAllInfoFieldAnnotations()));
            return union.stream().map(a -> a.getClass().getSimpleName()).collect(Collectors.toList());
        }

        /**
         * Annotation group names are simple names of all interfaces that are subtypes of {@ Annotation}.
         */
        public static List<String> getAllAnnotationGroupNames() {
            return ClassUtils.knownSubInterfaceSimpleNames(Annotation.class);
        }

        public List<InfoFieldAnnotation> createInfoFieldAnnotations() {
            final List<InfoFieldAnnotation> all = makeAllInfoFieldAnnotations();
            return filterAnnotations(all);
        }

        private static List<InfoFieldAnnotation> makeAllInfoFieldAnnotations() {
            return ClassUtils.makeInstancesOfSubclasses(InfoFieldAnnotation.class, Annotation.class.getPackage());
        }

        public List<GenotypeAnnotation> createGenotypeAnnotations() {
            final List<GenotypeAnnotation> all = makeAllGenotypeAnnotations();
            return filterAnnotations(all);
        }

        private static List<GenotypeAnnotation> makeAllGenotypeAnnotations() {
            return ClassUtils.makeInstancesOfSubclasses(GenotypeAnnotation.class, Annotation.class.getPackage());
        }

        /**
         * Returns a list of annotations that either:
         *  - belong to at least one of the requested annotation groups
         *  - belong to the set of requested annotations
         *
         *  - and are NOT listed for exclusion
         *
         *  The list is sorted by simple name of the class.
         */
        private <T extends VariantAnnotation> List<T> filterAnnotations(final List<T> all) {
            final SortedSet<T> annotations = new TreeSet<>(Comparator.comparing(t -> t.getClass().getSimpleName()));

            final Set<Class<?>> knownAnnotationGroups = ClassUtils.knownSubInterfaces(Annotation.class);

            for (final T t : all){
                if (!annotationsToExclude.contains(t.getClass().getSimpleName())) {
                    //if any group matches requested groups, it's in
                    @SuppressWarnings("unchecked")
                    final Set<Class<?>> annotationGroupsForT = ReflectionUtils.getAllSuperTypes(t.getClass(), sup -> sup.isInterface() && knownAnnotationGroups.contains(sup));
                    if (annotationGroupsForT.stream().anyMatch(iface -> annotationGroupsToUse.contains(iface.getSimpleName()))) {
                        annotations.add(t);
                    } else if (annotationsToUse.contains(t.getClass().getSimpleName())) {
                        annotations.add(t);
                    }
                }
            }

            return Collections.unmodifiableList(new ArrayList<>(annotations));
        }

    }
    protected static class VAExpression {

        public String fullName, fieldName;
        public FeatureInput<VariantContext> binding;
        public VCFInfoHeaderLine hInfo;

        public VAExpression(String fullExpression, List<FeatureInput<VariantContext>> dataSourceList){
            final int indexOfDot = fullExpression.lastIndexOf(".");
            if ( indexOfDot == -1 ) {
                throw new UserException.BadInput("The requested expression '"+fullExpression+"' is invalid, it should be in VCFFile.value format");
            }

            fullName = fullExpression;
            fieldName = fullExpression.substring(indexOfDot+1);

            final String bindingName = fullExpression.substring(0, indexOfDot);
            for ( final FeatureInput<VariantContext> ds : dataSourceList ) {
                if ( ds.getName().equals(bindingName) ) {
                    binding = ds;
                    break;
                }
            }
        }
    }

    protected List<VAExpression> getRequestedExpressions() { return expressions; }

    // select specific expressions to use
    public void addExpressions(Set<String> expressionsToUse, List<FeatureInput<VariantContext>> dataSources, boolean expressionAlleleConcordance) {//, Set<VCFHeaderLines>) {
        // set up the expressions
        for ( final String expression : expressionsToUse ) {
            expressions.add(new VAExpression(expression, dataSources));
        }
        this.expressionAlleleConcordance = expressionAlleleConcordance;
    }

    /**
     * Handles logic for expressions for variant contexts. Used to add annotations from one vcf file into the fields
     * of another if the variant contexts match sufficiently between the two files.
     *
     * @param vc  VariantContext to add annotations to
     * @param features  FeatureContext object containing extra VCF features to add to vc
     * @param ref  Reference context object corresponding to the region overlapping vc
     * @param attributes  running list of attributes into which to place new annotations
     */
    private void annotateExpressions(final VariantContext vc,
                                     final FeatureContext features,
                                     final ReferenceContext ref,
                                     final Map<String, Object> attributes){
        Utils.nonNull(vc);

        // each requested expression
        for ( final VAExpression expression : expressions ) {
            List<VariantContext> variantContexts = features.getValues(expression.binding, vc.getStart());

            if (!variantContexts.isEmpty()) {
                // get the expression's variant context
                VariantContext expressionVC = variantContexts.iterator().next();

                // special-case the ID field
                if (expression.fieldName.equals("ID")) {
                    if (expressionVC.hasID()) {
                        attributes.put(expression.fullName, expressionVC.getID());
                    }
                } else if (expression.fieldName.equals("ALT")) {
                    attributes.put(expression.fullName, expressionVC.getAlternateAllele(0).getDisplayString());
                } else if (expression.fieldName.equals("FILTER")) {
                    if (expressionVC.isFiltered()) {
                        attributes.put(expression.fullName, expressionVC.getFilters().toString().replace("[", "").replace("]", "").replace(" ", ""));
                    } else {
                        attributes.put(expression.fullName, "PASS");
                    }
                } else if (expressionVC.hasAttribute(expression.fieldName)) {
                    // find the info field
                    final VCFInfoHeaderLine hInfo = expression.hInfo;
                    if (hInfo == null) {
                        throw new UserException("Cannot annotate expression " + expression.fullName + " at " + ref.getInterval() + " for variant allele(s) " + vc.getAlleles() + ", missing header info");
                    }

                    //
                    // Add the info field annotations
                    //
                    final boolean useRefAndAltAlleles = VCFHeaderLineCount.R == hInfo.getCountType();
                    final boolean useAltAlleles = VCFHeaderLineCount.A == hInfo.getCountType();

                    // Annotation uses ref and/or alt alleles or enforce allele concordance
                    if ((useAltAlleles || useRefAndAltAlleles) || expressionAlleleConcordance) {

                        // remove brackets and spaces from expression value
                        final String cleanedExpressionValue = expressionVC.getAttribute(expression.fieldName,"").toString().replaceAll("[\\[\\]\\s]", "");

                        // get comma separated expression values
                        final ArrayList<String> expressionValuesList = new ArrayList<>(Arrays.asList(cleanedExpressionValue.split(",")));

                        boolean canAnnotate = false;
                        // get the minimum biallelics without genotypes

                        final List<VariantContext> minBiallelicVCs = getMinRepresentationBiallelics(vc);
                        final List<VariantContext> minBiallelicExprVCs = getMinRepresentationBiallelics(expressionVC);

                        // check concordance
                        final List<String> annotationValues = new ArrayList<>();
                        for (final VariantContext biallelicVC : minBiallelicVCs) {
                            // check that ref and alt alleles are the same
                            List<Allele> exprAlleles = biallelicVC.getAlleles();
                            boolean isAlleleConcordant = false;
                            int i = 0;
                            for (final VariantContext biallelicExprVC : minBiallelicExprVCs) {
                                List<Allele> alleles = biallelicExprVC.getAlleles();
                                // concordant
                                if (alleles.equals(exprAlleles)) {
                                    // get the value for the reference if needed.
                                    if (i == 0 && useRefAndAltAlleles) {
                                        annotationValues.add(expressionValuesList.get(i++));
                                    }
                                    // use annotation expression and add to vc
                                    annotationValues.add(expressionValuesList.get(i));
                                    isAlleleConcordant = true;
                                    canAnnotate = true;
                                    break;
                                }
                                i++;
                            }

                            // can not find allele match so set to annotation value to zero
                            if (!isAlleleConcordant) {
                                annotationValues.add("0");
                            }
                        }

                        // some allele matches so add the annotation values
                        if (canAnnotate) {
                            attributes.put(expression.fullName, annotationValues);
                        }
                    } else {
                        // use all of the expression values
                        attributes.put(expression.fullName, expressionVC.getAttribute(expression.fieldName));
                    }
                }
            }
        }
    }

    /**
     * Break the variant context into bialleles (reference and alternate alleles) and trim to a minimum representation
     *
     * @param vc variant context to annotate
     * @return list of biallelics trimmed to a minimum representation
     */
    private List<VariantContext> getMinRepresentationBiallelics(final VariantContext vc) {
        final List<VariantContext> minRepresentationBiallelicVCs = new ArrayList<>();
        if (vc.getNAlleles() > 2) {
            // TODO, this doesn't actually need to be done, we can simulate it at less cost
            for (int i = 1; i < vc.getNAlleles(); i++) {
                // Determining if the biallelic would have been considered a SNP
                if (! (vc.getReference().length() == 1 && vc.getAlternateAllele(i-1).length() == 1) ) {
                    minRepresentationBiallelicVCs.add(GATKVariantContextUtils.trimAlleles(
                            new VariantContextBuilder(vc)
                                    .alleles(Arrays.asList(vc.getReference(),vc.getAlternateAllele(i-1)))
                                    .attributes(removeIrrelevantAttributes(vc.getAttributes())).make(), true, true));
                } else {
                    minRepresentationBiallelicVCs.add(new VariantContextBuilder(vc)
                            .alleles(Arrays.asList(vc.getReference(),vc.getAlternateAllele(i-1)))
                            .attributes(removeIrrelevantAttributes(vc.getAttributes())).make());
                }
            }
        } else {
            minRepresentationBiallelicVCs.add(vc);
        }

        return minRepresentationBiallelicVCs;
    }

    private Map<String, Object> removeIrrelevantAttributes(Map<String, Object> attributes) {
        // since the VC has been subset, remove the invalid attributes
        Map<String, Object> ret = new HashMap<>(attributes);
        for ( final String key : attributes.keySet() ) {
            if ( !(key.equals(VCFConstants.ALLELE_COUNT_KEY) || key.equals(VCFConstants.ALLELE_FREQUENCY_KEY) || key.equals(VCFConstants.ALLELE_NUMBER_KEY)) ) {
                ret.remove(key);
            }
        }

        return ret;
    }

}
