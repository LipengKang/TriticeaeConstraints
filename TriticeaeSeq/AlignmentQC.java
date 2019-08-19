/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import lipengKang.analysis.KStringUtils;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author kanglipeng
 */
//Triticeae PacBio Seq data QC
public class AlignmentQC {

    int sum = 0;

    public AlignmentQC() {
//this.orthOverlapRate();
         this.geneCoverRate1(); //test maf cover rate of gene
//this.geneCoverRate2();
//this.geneCoverRate3();
       // this.geneCoverRate4();
    }

    public void orthOverlapRate() {
        String bed = "/data1/home/lipeng/software/mummer-4.0.0beta2/wheatA.barley.coords";
        String maf = "/data1/home/lipeng/result/TriticeaeCM/align/A/wheatA.tu-ctg.split2.maf";
        // String bed="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.coords";
        //String maf="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.maf";
        String temp = null;
        String tem[] = null;
        BufferedReader br1 = IOUtils.getTextReader(bed);
        BufferedReader br2 = IOUtils.getTextReader(maf);
        int[][] nucmerIndexArr = new int[7][800000000];

        long nucmerCoverScore = 0;

        try {
            temp = br1.readLine();
            temp = br1.readLine();
            temp = br1.readLine();
            temp = br1.readLine();
            while ((temp = br1.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                for (int q = Integer.parseInt(tem[0]); q <= Integer.parseInt(tem[1]); q++) {
                    nucmerIndexArr[Integer.valueOf(tem[7].substring(10)) - 1][q] = 1;

                }
            }
            br1.close();
            //  int geneCoverScore=0;
            int[] overlapScore = new int[7];
            long sumOverlapCover = 0;
            long mafCoverScore = 0;
            for (int i = 0; i < 7; i++) {
                for (int j = 0; j < nucmerIndexArr[i].length; j++) {
                    nucmerCoverScore = nucmerCoverScore + nucmerIndexArr[i][j];
                }
            }

            while ((temp = br2.readLine()) != null) {
                if (!temp.startsWith("s wheatA")) {
                    continue;
                }
                List<String> fList = KStringUtils.fastSplit(temp);
                List<String> fListNew = new ArrayList<>();
                for (int i = 0; i < fList.size(); i++) {
                    if (fList.get(i) != null && !fList.get(i).equals("")) {
                        fListNew.add(fList.get(i));
                    }
                }
                tem = fListNew.toArray(new String[fListNew.size()]);

                for (int p = Integer.parseInt(tem[2]) + 1; p <= Integer.parseInt(tem[2]) + Integer.parseInt(tem[3]); p++) {
                    mafCoverScore++;
                    if (nucmerIndexArr[Integer.valueOf(tem[1].substring(10)) - 1][p] != 1) {
                        continue;
                    }
// geneCoverScore= geneCoverScore+geneIndexArr[Integer.valueOf(tem[1].substring(2))-1][p];
                    overlapScore[Integer.valueOf(tem[1].substring(10)) - 1] = overlapScore[Integer.valueOf(tem[1].substring(10)) - 1] + 1;
                }
            }
            br2.close();

            for (int i : overlapScore) {
                System.out.println("maf covers " + i + " bp chromosome genic region");
                sumOverlapCover = sumOverlapCover + i;
            }
            System.out.println("maf covers " + sumOverlapCover + " bp overlap region");
            //  System.out.println("maf covers " +geneCoverScore+" bp genic region");
            //  System.out.println("maf covers " +geneCoverScore+" bp genic region");
            System.out.println("nucmer covers " + nucmerCoverScore + " bp genomic region");
            System.out.println("maf covers " + mafCoverScore + " bp genomic region");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    //-----------------------------test maf cover rate of gene  --------------------------
    public void geneCoverRate1() {
        String bed = "/data1/home/lipeng/database/genome/wheat/annotation/iwgsc_refseqv1.1_genes_2017July06/CNEfinderGene.gff3";
        String maf = "/data1/home/lipeng/result/TriticeaeCM/pairWiseAlign/T_aestivum-A_tauschii.w3.sing.maf";
        String temp = null;
        String tem[] = null;
        BufferedReader br1 = IOUtils.getTextReader(bed);
        BufferedReader br2 = IOUtils.getTextReader(maf);
        int[][] geneIndexArr = new int[7][800000000];
int geneticLength=0;
        try {
            br1.readLine();
            while ((temp = br1.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if(!tem[1].equals("chr1A"))continue;
                for (int q = Integer.parseInt(tem[2]); q <= Integer.parseInt(tem[3]); q++) {
                    geneIndexArr[Integer.parseInt(String.valueOf(tem[1].charAt(3))) - 1][q] = 1;
                }
                geneticLength =geneticLength+Integer.parseInt(tem[3])-Integer.parseInt(tem[2]);
            }
            br1.close();
            //  int geneCoverScore=0;
            int[] geneCoverScore = new int[7];
            long geneSumCover = 0;
            long coverScore = 0;

            while ((temp = br2.readLine()) != null) {
                if (!temp.startsWith("s Tr")) {
                    continue;
                }
                List<String> fList = KStringUtils.fastSplit(temp);
                List<String> fListNew = new ArrayList<>();
                for (int i = 0; i < fList.size(); i++) {
                    if (fList.get(i) != null && !fList.get(i).equals("")) {
                        fListNew.add(fList.get(i));
                    }
                }
                tem = fListNew.toArray(new String[fListNew.size()]);

                for (int p = Integer.parseInt(tem[2]) + 1; p <= Integer.parseInt(tem[2]) + Integer.parseInt(tem[3]); p++) {
                    coverScore++;
                    if (geneIndexArr[Integer.valueOf(String.valueOf(tem[1].charAt(9)))- 1][p] != 1) {
                        continue;
                    }
// geneCoverScore= geneCoverScore+geneIndexArr[Integer.valueOf(tem[1].substring(2))-1][p];
                    geneCoverScore[Integer.valueOf(String.valueOf(tem[1].charAt(9))) - 1] = geneCoverScore[Integer.valueOf(String.valueOf(tem[1].charAt(9))) - 1] + 1;
                }
            }
            br2.close();

            for (int i : geneCoverScore) {
                //System.out.println("maf covers " + i + " bp chromosome genic region");
                geneSumCover = geneSumCover + i;
            }
            System.out.println("maf covers " + geneSumCover + " bp genic region");
            System.out.println("maf covers " +100*geneSumCover/geneticLength+" % genic region");
            System.out.println("maf covers " + coverScore + " bp genomic region");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void geneCoverRate2() {
        String bed = "/data1/home/lipeng/result/PacBioAnalysis/maptest/10x.sam";
        // String bed="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.sam";
        String temp = null;
        String tem[] = null;
        BufferedReader br1 = IOUtils.getTextReader(bed);
        int[][] genomeIndexArr = new int[7][760000000];

        try {
            while ((temp = br1.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (!tem[2].startsWith("Tu")) {
                    continue;
                }
                for (int q = Integer.parseInt(tem[3]); q <= Integer.parseInt(tem[3]) + tem[9].length() - 1; q++) {
                    genomeIndexArr[Integer.valueOf(tem[2].substring(2)) - 1][q] = 1;
                }
            }
            br1.close();
            long coverLength = 0;
            for (int i = 0; i < 7; i++) {
                for (int j = 0; j < genomeIndexArr[i].length; j++) {
                    coverLength = coverLength + genomeIndexArr[i][j];
                }
            }

            System.out.println("maf covers " + coverLength + " bp genomic region");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void geneCoverRate3() {

        String maf = "/data1/home/lipeng/result/PacBioAnalysis/homo6/RAC.postmask.maf";
        String temp = null;
        String tem[] = null;

        BufferedReader br1 = IOUtils.getTextReader(maf);

        int[] geneIndexArr = new int[258957000];
        String bed = "/data1/home/lipeng/software/CNEFinder/Experiments/Files/hg38_genes";

        BufferedReader br2 = IOUtils.getTextReader(bed);

        try {
            while ((temp = br2.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (!tem[1].equals("1")) {
                    continue;
                }
                for (int q = Integer.parseInt(tem[2]); q <= Integer.parseInt(tem[3]); q++) {
                    geneIndexArr[q] = 1;
                }
            }
            br2.close();
            int geneCoverBps = 0;
            while ((temp = br1.readLine()) != null) {
                if (!temp.startsWith("s chr")) {
                    continue;
                }
                List<String> fList = KStringUtils.fastSplit(temp);
                List<String> fListNew = new ArrayList<>();
                for (int i = 0; i < fList.size(); i++) {
                    if (fList.get(i) != null && !fList.get(i).equals("")) {
                        fListNew.add(fList.get(i));
                    }
                }
                tem = fListNew.toArray(new String[fListNew.size()]);

                for (int p = Integer.parseInt(tem[2]) + 1; p <= Integer.parseInt(tem[2]) + Integer.parseInt(tem[3]); p++) {
                    if (geneIndexArr[p] != 1) {
                        continue;
                    }
// geneCoverScore= geneCoverScore+geneIndexArr[Integer.valueOf(tem[1].substring(2))-1][p];
                    geneCoverBps = geneCoverBps + 1;
//geneIndexArr[p]=1;
                }
            }
            br1.close();
            int geneBps = 0;
            for (int i : geneIndexArr) {
                geneBps = geneBps + i;
            }
//   int coverBps=0;
//    for (int i: geneIndexArr){
//    coverBps=coverBps+i;
//    }
//    System.out.println("maf covers " +coverBps+" bp chr1");
            System.out.println("maf covers " + geneCoverBps + " bp chr1 gene");
            System.out.println("chr1 owns " + geneBps + " bp gene region");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void geneCoverRate4() {
        String bed = "/data1/home/lipeng/result/PacBioAnalysis/homo6/aln.sam";
        // String bed="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.sam";
        String temp = null;
        String tem[] = null;
        BufferedReader br1 = IOUtils.getTextReader(bed);
        int[] genomeIndexArr = new int[258956422];

        try {
            while ((temp = br1.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);

                for (int q = Integer.parseInt(tem[7]); q <= Integer.parseInt(tem[8]); q++) {
                    genomeIndexArr[q] = 1;
                }
            }
            br1.close();
            long coverLength = 0;

            for (int j = 0; j < genomeIndexArr.length; j++) {
                coverLength = coverLength + genomeIndexArr[j];
            }

            System.out.println("paf covers " + coverLength + " bp genomic region");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        new AlignmentQC();
    }
}
