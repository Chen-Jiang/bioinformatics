import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class FindQualifiedVariants {

    ArrayList<VariantContext> variantList = new ArrayList<>();
    ArrayList<VariantContext> variantAndInsideExonList = new ArrayList<>();
    ReadExonFile readEx = new ReadExonFile();
    CalculateVariantsPossibility calculateVariantsPossibility = new CalculateVariantsPossibility();

    public void Read() {

        // ============step1: read the family vcf file to extract useful information============

        System.out.println("For Sickle Cell Anaemia:");
        // reader1 just for test.....
        VCFFileReader reader1 = new VCFFileReader(new File("/Users/Shawn/Desktop/t1.vcf.gz"));
        VCFFileReader reader = new VCFFileReader(new File("/Users/Shawn/Downloads/merged.vcf.gz"));
        CloseableIterator<VariantContext> text = reader.iterator();

        int num = 0;

        while (text.hasNext()) {
            VariantContext context = text.next();

            // alleles is a list containing both reference gene and alt gene
            List<Allele> alleles = context.getAlleles();

            // get this record's reference allele, generally ends with "*"
            String reference = context.getReference().toString();
            if (reference.endsWith("*")) {
                reference = reference.substring(0, reference.length()-1);
            }

            // get this record's alt allele
            List<Allele> alt = context.getAlternateAlleles();

            int count = 0;

            // get all the seven family members' information
            // the family members follow the order: father, mother, daughter1, d2, d3, son1, s2
            for (int i = 0; i < 7; i++) {
                String name = context.getGenotype(i).getSampleName();
                String genotype = context.getGenotype(i).getGenotypeString();

                // look for qualified variants that match the inheritance pattern

                // parents, daughter1 and son1 has one reference allele and one slt allele
                if (i == 0 || i == 1 || i == 2 || i == 5) {
                    if (!(genotype.substring(0, 1).equals(genotype.substring(genotype.length()-1)))) {
                        if (genotype.contains(reference)) {
                            for (int j = 0; j< alt.size(); j++) {
                                if (genotype.contains(alt.get(j).toString())) {
                                    count++;
                                }
                            }

                        }
                    }
                }
                // daughter2 and son2 has two alt alleles
                if (i == 3 || i == 6) {
                    if ((genotype.substring(0, 1).equals((genotype.substring(genotype.length()-1)))) && !(genotype.contains(reference))) {
                        count++;
                    }

                }
                //daughter3 has two reference alleles
                if (i == 4) {
                    if ((genotype.substring(0, 1).equals((genotype.substring(genotype.length()-1)))) && genotype.contains(reference)) {
                        count++;
                    }
                }
                if (count != i+1) {
                    break;
                }
            }
            // put all the qualified variants into a new list
            if (count == 7) {
                variantList.add(context);
            }
        }
        System.out.println("size: " + variantList.size());



        // ====step2: find variants which locates inside exonic regions among all the variants in the variantList=====

        for (VariantContext item : variantList) {
            if (readEx.checkExonRegion(item)) {
                variantAndInsideExonList.add(item);
            }
        }
        System.out.println("new size: " + variantAndInsideExonList.size());


        // ============step4: calculate the possibility of variants according to this file:============
        // ============" ALL.merged.phase1_release_v3.20101123.snps_indels_svs.vcf "============

        for (VariantContext item: variantAndInsideExonList) {
            double variantPossibility = 0;
            variantPossibility = calculateVariantsPossibility.calculateVariantPossibility(item);
            System.out.println(item.getContig() + "--" + item.getStart());
            System.out.println(item.getAlleles());
            for (int i = 0; i < 7; i++) {
                System.out.println(item.getGenotype(i).getSampleName() + ": " + item.getGenotype(i).getGenotypeString());
            }
            System.out.println("possibility: " + variantPossibility);
            System.out.println("===========");
        }
        text.close();
        reader.close();
    }

    public static void main(String[] args) {
        new FindQualifiedVariants().Read();
    }
}
