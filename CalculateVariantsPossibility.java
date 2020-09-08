import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CalculateVariantsPossibility {


    /**
     * This is to read all people's variants file: "ALL.merged.phase1_release_v3.20101123.snps_indels_svs.vcf.gz"
     * @return a list containing useful information
     */

    public ArrayList<Map<String, String>> read() {

        ArrayList<Map<String, String>> variantsList = new ArrayList<>();

        // reader1 just for test.....
        VCFFileReader reader1 = new VCFFileReader(new File("/Users/Shawn/Desktop/t2.vcf.gz"));
        VCFFileReader reader = new VCFFileReader(new File("/Users/Shawn/Downloads/ALL.merged.phase1_release_v3.20101123.snps_indels_svs.vcf.gz"));
        CloseableIterator<VariantContext> text = reader.iterator();

        while (text.hasNext()) {
            Map<String, String> variant = new HashMap<>();

            List<Allele> alleles = new ArrayList<>();
            VariantContext context = text.next();
            String chromo = context.getContig();
            String start = context.getStart()+"";
            String end = context.getEnd()+"";
            String reference = context.getReference().toString();
            String alt = context.getAlternateAlleles().toString();

            // Alternate Allele Count
            String AC = context.getAttribute("AC").toString();
            // Total Allele Count
            String AN = context.getAttribute("AN").toString();

            variant.put("chromo", chromo);
            variant.put("start",start);
            variant.put("end", end);
            variant.put("reference", reference);
            variant.put("alt",alt);
            variant.put("AC", AC);
            variant.put("AN",AN);

            variantsList.add(variant);
        }
        return variantsList;
    }

    /**
     * Calculate the variants possibility according to AC and AN
     * @param variantContext one single variant record
     * @return the possibility
     */

    public double calculateVariantPossibility(VariantContext variantContext) {
        ArrayList<Map<String, String>> variantsList = read();
        double variantPossibility = 0;
        ArrayList<Map<String, String>> sameChromo = new ArrayList<>();
        int num = 0;

        // get the information: which is the current chromosomal
        String chromo = variantContext.getContig();

        // find the same chromosomal from the list read from sample population VCF data
        for (Map<String, String> item: variantsList) {
            if (item.get("chromo").equals(chromo)) {
                sameChromo.add(item);
            }
        }

        // find all the same chromosomal records to see if there is a same position
        // if has the same position, calculate the AC/AN
        for (Map<String, String> item: sameChromo) {
            if (Integer.parseInt(item.get("start")) == (variantContext.getStart())) {
                num++;
                double AC = Double.parseDouble(item.get("AC"));
                double AN = Double.parseDouble(item.get("AN"));
                variantPossibility = AC/AN;
            }
        }

        // if no same position found in the list, add the family member to the list and calculate
        if (num == 0) {
            int count = 0;
            double AN = 2184.0;
            for (int j = 0; j<7; j++) {
                List<Allele> alt = variantContext.getAlternateAlleles();
                for (int k = 0; k<alt.size(); k++) {
                    if (variantContext.getGenotype(j).getGenotypeString().contains(alt.get(k).toString())) {
                        count++;
                    }
                }
            }
            variantPossibility = count/(AN+count);
        }
        return variantPossibility;
    }


    public static void main(String[] args) {
        new CalculateVariantsPossibility().read();
    }
}
