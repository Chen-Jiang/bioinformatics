# bioinformatics
This is an assignment of one of my course, Bioinformatics. In the assignment, we use a variant call analysis method to find the possible genes which can cause family inherited diseases.  Considering the huge amount of human genes, we use programming to help us find the genes more quickly.
All the steps included in the files are below:
1. Read the exon file and put all the exon regions to a Map;
2. Found all the positions from family VCF file which fit the family disease inheritance pattern, and use the exon Map as a filter, to remove all the variants which are not located in exon regions;
3. Read the population VCF file which contains 1092 people’s variants, and compare the variants’ chromosome and position information we found from family VCF to see if the variant is also found in the 1092 persons’ DNA. If not found, add this record;
4. Explore the prevalence of each disease;
