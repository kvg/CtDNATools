package org.broadinstitute.hellbender.utils.codecs.GENCODE;

import java.util.ArrayList;

/**
 * A Gencode GTF Feature representing an untranslated region.
 *
 * A GTF Feature represents one row of a GTF File.
 * The specification of a GTF file is defined here:
 * http://mblab.wustl.edu/GTF22.html
 *
 * Created by jonn on 7/25/17.
 */
final public class GencodeGtfUTRFeature extends GencodeGtfFeature {

    private GencodeGtfUTRFeature(String[] gtfFields) {
        super(gtfFields);
    }

    public static GencodeGtfFeature create(String[] gtfFields) {
        return new GencodeGtfUTRFeature(gtfFields);
    }

    private GencodeGtfUTRFeature(long featureOrderNumber,
                                 String chromosomeName,
                                 AnnotationSource annotationSource,
                                 FeatureType featureType,
                                 int genomicStartLocation,
                                 int genomicEndLocation,
                                 GenomicStrand genomicStrand,
                                 GenomicPhase genomicPhase,
                                 String geneId,
                                 String transcriptId,
                                 GeneTranscriptType geneType,
                                 GeneTranscriptStatus geneStatus,
                                 String geneName,
                                 GeneTranscriptType transcriptType,
                                 GeneTranscriptStatus transcriptStatus,
                                 String transcriptName,
                                 int exonNumber,
                                 String exonId,
                                 LocusLevel locusLevel,
                                 ArrayList<OptionalField<?>> optionalFields,
                                 String anonymousOptionalFields) {

        super(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
    }

    public static GencodeGtfFeature create(long featureOrderNumber,
                                           String chromosomeName,
                                           AnnotationSource annotationSource,
                                           FeatureType featureType,
                                           int genomicStartLocation,
                                           int genomicEndLocation,
                                           GenomicStrand genomicStrand,
                                           GenomicPhase genomicPhase,
                                           String geneId,
                                           String transcriptId,
                                           GeneTranscriptType geneType,
                                           GeneTranscriptStatus geneStatus,
                                           String geneName,
                                           GeneTranscriptType transcriptType,
                                           GeneTranscriptStatus transcriptStatus,
                                           String transcriptName,
                                           int exonNumber,
                                           String exonId,
                                           LocusLevel locusLevel,
                                           ArrayList<OptionalField<?>> optionalFields,
                                           String anonymousOptionalFields) {

        return new GencodeGtfUTRFeature(featureOrderNumber, chromosomeName, annotationSource, featureType, genomicStartLocation, genomicEndLocation, genomicStrand, genomicPhase, geneId, transcriptId, geneType, geneStatus, geneName, transcriptType, transcriptStatus, transcriptName, exonNumber, exonId, locusLevel, optionalFields, anonymousOptionalFields);
    }

}
