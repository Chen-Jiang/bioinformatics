import htsjdk.variant.variantcontext.VariantContext;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class ReadExonFile {


    /**
     *  This is to read the exon file, the file name is : wgEncodeGencodeBasicV17.txt
     */
    public ArrayList<Map<String, String>> read() {

        ArrayList<Map<String, String>> exonList = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(new File("/Users/Shawn/Desktop/VCF/wgEncodeGencodeBasicV17.txt")));){

            String line = null;
            while ((line = reader.readLine()) != null) {
                Map<String, String> exon = new HashMap<>();

                // this txt file is a tab-delimited format
                String[] va = line.split("\t");

                // get chromosomal, number of exons, start position and end position fo each exon region
                String chromo = va[2].substring(3);
                String numExon = va[8];
                String exonStarts = va[9];
                String exonEnds = va[10];

                exon.put("chromo",chromo);
                exon.put("num", numExon);
                exon.put("starts", exonStarts);
                exon.put("ends", exonEnds);

                exonList.add(exon);
            }

        } catch (IOException exception) {
            exception.printStackTrace();
        }
        return exonList;
    }

    public boolean checkExonRegion(VariantContext variantContext) {

        Boolean isInsideExonRegon = false;

        ArrayList<Map<String, String>> exonList = read();

        for (Map<String, String> item: exonList) {
            String starts = item.get("starts");
            String ends = item.get("ends");
            int numExon = Integer.parseInt(item.get("num"));

            String[] exonStarts = starts.split(",");
            String[] exonEnds = ends.split(",");

            // to see if in the same chromosomal
            if (variantContext.getContig().equals(item.get("chromo"))) {
                int position = variantContext.getStart();
                // to see if included one of the exon regions
                for (int i = 0; i<numExon; i++) {
                    if (position >= Integer.parseInt(exonStarts[i]) && position <= Integer.parseInt(exonEnds[i])) {
                        isInsideExonRegon = true;
                    }
                }
            }
        }
        return isInsideExonRegon;
    }

    public static void main(String[] args) {
        new ReadExonFile().read();
    }
}
