/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import com.koloboke.collect.set.hash.HashIntSet;
import static com.koloboke.collect.set.hash.HashIntSets.newMutableSet;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import static java.util.Collections.list;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import java.util.List;
import java.util.Set;
import lipengKang.analysis.KStringUtils;
import lipengKang.analysis.range;

import utils.IOUtils;
import utils.PStringUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author kanglipeng
 */
//wheatGeneQual aims to call wheat gene basic situation(such as length distribution, etc) 
//geneReadsCounter aims to count PacBio reads mapping quality in gene
public class wheatGeneQual {

    public wheatGeneQual() {
        //  this.countXGene();
       // this.callGeneLength();//this function input is standard IWSC gff3 ,not gff3 modified by fei
        // this.indexGene();     //index genes
        // this.parseGerpelem();
        // this.callMappedGenePos();
        // this.callGeneMappedDepth(outFile);           //need update .......
        //this.callChrElementPos();  //create position lists for UTRs,exon,gene from gff3
        //this.extratLongestTranscript();
        //this.callgeneticElementLength();
        // this.callMappedDepth(depthFile, subGenome); // need update.......
        //this.intrageneGERP();          // (require callChrElementPos) GERP distribution in intragenetic regions
        //this.callGERPConstraintNum(subGenome);    //call GERP >? regions
        //this.GERPWindowDistribution();            //GERP ?bp window distribution of chromosomes
         //this.figMAF();
         //this.sFSofDeletrious();
     //   this.deleteriousstat();
       
    }
    String gff3Dir = null;
    String subGenome = null;
    String outFile = null;
    String samFile = null;
    String depthFile = null;
    String temp = null;
    String[] tem = null;
    TIntArrayList mappedPos = new TIntArrayList();
    HashMap<Integer, HashIntSet> mappedReadsPosMap = new HashMap<>();
    HashMap<Integer, String> chrGeneIndexMap = new HashMap<>();
    HashMap<Integer, HashMap<Integer, String>> allGeneIndexMap = new HashMap<>();
    ArrayList<Integer> genePosList = new ArrayList();
    HashMap<Integer, HashIntSet> allGenePosMap = new HashMap<>();
    HashMap<String, Integer> mappedDepthMap = new HashMap<>();
    HashMap<Integer, HashMap<String, ArrayList<Integer>>> ElementPosMap = new HashMap();
    ArrayList<String> readsNameList = new ArrayList<>();

//--------------------------------------call gene or intrageneic elements length or numbers from gff3------------------------
    //used to following hist plot gene length
    private void callGeneLength() {
        // String gff3Dir = "/data1/home/lipeng/database/genome/urartu/WheatTu.annotation/WheatTu.gene.gff";
        String gff3Dir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/chrgene.bed";
        BufferedReader br;
        //   BufferedWriter bw;
        if (gff3Dir.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(gff3Dir);
            //     bw = YaoIOUtils.getTextGzipWriter(outFile);
        } else {
            br = YaoIOUtils.getTextReader(gff3Dir);
            //       bw = YaoIOUtils.getTextWriter(outFile);
        }
        int StartPos = 0;
        int EndPos = 0;
        long Length = 0;
        long geneLength = 0;
//        long exonLength = 0;
//        long fiveUTRLength = 0;
//        long threeUTRLength = 0;
        //int Num = 0;
        //    ArrayList<Integer> geneLengthList = new ArrayList<>();
        try {

            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("Tu")) {
                    List<String> tList = PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);

                    if (tem[7].equals("gene") && tem[0].length() == 3) {
                        StartPos = Integer.parseInt(tem[1]);
                        EndPos = Integer.parseInt(tem[2]);
                        Length = EndPos - StartPos;
                        //  Num++;
                        //  geneLengthList.add(geneLength);
                        geneLength = geneLength + Length;
                    }
                    /* if (tem[2].equals("exon") && tem[0].endsWith(subGenome)) {
                        StartPos = Integer.parseInt(tem[3]);
                        EndPos = Integer.parseInt(tem[4]);
                        Length = EndPos - StartPos;
                        exonLength = exonLength + Length;
                    }
                    if (tem[2].equals("five_prime_UTR") && tem[0].endsWith(subGenome)) {
                        StartPos = Integer.parseInt(tem[3]);
                        EndPos = Integer.parseInt(tem[4]);
                        Length = EndPos - StartPos;
                        fiveUTRLength = fiveUTRLength + Length;
                    }
                    if (tem[2].equals("three_prime_UTR") && tem[0].endsWith(subGenome)) {
                        StartPos = Integer.parseInt(tem[3]);
                        EndPos = Integer.parseInt(tem[4]);
                        Length = EndPos - StartPos;
                        threeUTRLength = threeUTRLength + Length;
                    }*/

                }
            }
            System.out.println("urartu genes sum length : " + geneLength);
            /*System.out.println("wheat " + subGenome + "all exon length : " + exonLength);
            System.out.println("wheat " + subGenome + "all threeUTR length : " + threeUTRLength);
            System.out.println("wheat " + subGenome + "all fiveUTR length : " + fiveUTRLength);*/
            //System.out.print("There are" + " " + geneNum + " " + "genes in wheat" + subGenome + " genome");
            /* Collections.sort(geneLengthList);
            for (int i = 0; i < geneLengthList.size(); i++) {

                bw.write(String.valueOf(geneLengthList.get(i)));
                bw.newLine();
            }*/

            br.close();
            // bw.flush();
            // bw.close();

        } catch (Exception e) {
            System.out.println("Error in calling genes length of " + "urartu genome");
            e.printStackTrace();
        }
    }

    public void countXGene() {
        String pepDir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_LC_20170706_pep.fasta";
        BufferedReader br;
        if (pepDir.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(pepDir);
        } else {
            br = YaoIOUtils.getTextReader(pepDir);
        }
        try {

            int xgeneNum = 0;

            while ((temp = br.readLine()) != null) {
//                while(temp.startsWith(">")){
//                   temp = br.readLine();
//                  while(!temp.contains("X")) {
//                        temp = br.readLine();
//                        if(temp==null){break;}
//                        if(temp.startsWith(">")){
//                            break;
//                        }
//                      
//                    }
//                  if(temp==null){break;}
//                   if(temp.contains("X")){
//                   xgeneNum++;}
//                }
                if (temp.startsWith(">")) {
                    xgeneNum++;
                }

            }

            System.out.print(xgeneNum);
            br.close();

        } catch (Exception e) {
            System.out.println("Error in calling genes length of " + "urartu genome");
            e.printStackTrace();

        }
    }
    //--------------------------indexGene positions--------------------------
    //create index of all gene Pos in GFF3

    private HashMap<Integer, HashIntSet> indexGene() {
        String gff3Dir = "/data1/home/lipeng/database/genome/urartu/WheatTu.annotation/WheatTu.gene.gff";
        //  String gff3Dir ="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.gff3";
        //String gff3Dir="/data1/home/lipeng/database/genome/wheat/annotation/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gff3";
        BufferedReader br;
        if (gff3Dir.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(gff3Dir);
        } else {
            br = YaoIOUtils.getTextReader(gff3Dir);
        }
        int geneStartPos = 0;
        int geneEndPos = 0;
//        String geneInfo = null;

        HashIntSet[] genePosSetArr = new HashIntSet[7];
        for (int i = 0; i <= 6; i++) {
            genePosSetArr[i] = newMutableSet();
        }
//         HashMap[] geneIndexArr = new HashMap[7];

        try {

            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("TuUngrouped")) {
                    continue;
                }
// if (temp.startsWith("chr1A")) {

                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
//                    StringBuilder chr = new StringBuilder();
//                    chr.append(String.valueOf(i)).append("A");
                if (tem[2].equals("gene")) {
                    int chr = Integer.valueOf(String.valueOf(tem[0].charAt(2)));
//                    if(geneStartPos>=1000){
//                    geneStartPos = Integer.parseInt(tem[3])-10000;
//                    geneEndPos = Integer.parseInt(tem[4])+10000;}else{
//                    geneStartPos = 0;
//                     geneEndPos = Integer.parseInt(tem[4])+10000;
//                    }

                    geneStartPos = Integer.parseInt(tem[3]);
                    geneEndPos = Integer.parseInt(tem[4]);
                    // geneInfo = KStringUtils.fastSplitSemicolon(tem[8]).get(0).substring(3);
                    // hash all gene Position in a chr ------>gene name
                    for (int j = geneStartPos; j <= geneEndPos; j++) {
                        //  geneIndexArr[chr].put(j, geneInfo);
                        genePosSetArr[chr - 1].add(j);
                        // genePosList.add(j);
                    }
                }
            }

            br.close();
            for (int i = 1; i < 8; i++) {
                allGenePosMap.put(i, genePosSetArr[i - 1]);
                //  allGeneIndexMap.put(i, geneIndexArr[i]);
            }
//                genePosList = new ArrayList<Integer>();
//                chrGeneIndexMap = new HashMap<Integer, String>();
            System.out.println("indexing genes in urartu genome. ok!");
        } catch (Exception e) {
            System.out.println("Error in indexing genes in urartu genome");
            e.printStackTrace();
        }
        return allGenePosMap;
    }
//-----------------------------parse minimap2->.PAF ------------------------------------------
    //collect mapped contigs-gene pos 

    public void callMappedGenePos() {
        String pafFile = "/data1/home/lipeng/result/PacBioAnalysis/urartuQC1/readschain.kbmap";
        //   String pafFile ="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.paf";
        BufferedReader br;
        if (pafFile.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(pafFile);
        } else {
            br = YaoIOUtils.getTextReader(pafFile);
        }
        int mappedBlockStart = 0;
        int mappedBlockEnd = 0;
        long contigsSumLength = 0;
        HashIntSet[] mappedPosArr = new HashIntSet[7];
        for (int i = 0; i <= 6; i++) {
            mappedPosArr[i] = newMutableSet();
        }
        try {
            while ((temp = br.readLine()) != null) {

                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (tem[5].startsWith("TuU")) {
                    continue;
                }

                mappedBlockStart = Integer.parseInt(tem[8]);
                mappedBlockEnd = Integer.parseInt(tem[9]);
                contigsSumLength = contigsSumLength + mappedBlockEnd - mappedBlockStart;
                int chr = Integer.valueOf(String.valueOf(tem[5].charAt(2)));
                for (int j = mappedBlockStart; j <= mappedBlockEnd; j++) {
                    mappedPosArr[chr - 1].add(j);
//                    mappedPosList.add(j);
                }

            }
            System.out.println("indexing wtdbg2-contigs. ok!");
            br.close();
            //create mapped contigs pos lists for each chr 
            for (int i = 1; i <= 7; i++) {

                mappedReadsPosMap.put(i, mappedPosArr[i - 1]);
            }

            System.out.println(" wtdbg2-contigs index map constructing. ok!");

            //start find intersection between contigs and reference
            HashMap<Integer, HashIntSet> allGenePosMap = this.allGenePosMap;
            for (int j = 1; j <= 7; j++) {
//               ArrayList <Integer>interSectionPos = new ArrayList<>();
//               interSectionPos= range.getCommonElements(mappedPosList,allGenePosMap.get(j));
                HashIntSet chrGenePosSet = allGenePosMap.get(j);
                HashIntSet chrmappedPosSet = mappedReadsPosMap.get(j);
                chrmappedPosSet.retainAll(chrGenePosSet);

                System.out.println("Urartu chr " + j + " was covered " + chrmappedPosSet.size() + " bp by wtdbg2-contigs");
            }
            System.out.println("paf contigs cover " + contigsSumLength + " bp");

        } catch (Exception e) {
            System.out.println("Error in calling mapped reads position ");
            e.printStackTrace();
        }

    }

//    public void callGeneMappedDepth(String outFile) {
//        //// allGeneIndexMap:chrNum----->( gene Positions in a chr ------>gene name)
//        //mappedReadsPosMap:chrNum-gene numeric number----->mappedReadsPos
//        //allGenePosMap:chrNum------>gene Positions in a chr
//
//        HashMap<Integer, TIntArrayList> mappedReadsPosMap = this.mappedReadsPosMap;
//        HashMap<Integer, ArrayList<Integer>> allGenePosMap = this.allGenePosMap;
//        ArrayList<String> mappedReadsNameList = new ArrayList<>();
//        ArrayList<String> mappedGeneNameList = new ArrayList();
//        Set<String> mappedGeneNameSet = new HashSet();
//        ArrayList<Integer> mappedDepthList = new ArrayList<>();
//        for (int i = 0; i <= 7; i++) {
//            ArrayList<Integer> chrGenePosList = new ArrayList();
//            chrGenePosList = allGenePosMap.get(i);
//            Collections.sort(chrGenePosList);
//            int mappedDepth = 0;
//            int mappedBaseNum = 0;
//            for (int j = 0; j < mappedReadsNameList.size(); j++) {
//                if (mappedReadsNameList.get(j).startsWith(String.valueOf(i))) {
//                    for (int k = mappedReadsPosMap.get(mappedReadsNameList.get(j)).get(0); k <= mappedReadsPosMap.get(mappedReadsNameList.get(j)).get(1); k++) {
//                        if (Collections.binarySearch(chrGenePosList, k) >= 0) {
//                            mappedBaseNum++;
//                            if (mappedBaseNum >= 30) {
//                                mappedGeneNameSet.add(allGeneIndexMap.get(i).get(k));
//                            }
//                        }
//                        ArrayList<String> perReadsMappedGeneList = new ArrayList(mappedGeneNameSet);
//                        mappedGeneNameSet = new HashSet();
//                        for (int l = 0; l < perReadsMappedGeneList.size(); l++) {
//                            mappedGeneNameList.add(perReadsMappedGeneList.get(l));
//                        }
//                        mappedBaseNum = 0;
//
//                        //--------
//                        //      ---------   each seed match >=15bp is defined to mapped
//                    }
//
//                }
//            }
//        }
//        Collections.sort(mappedGeneNameList);
//        for (String gene : mappedGeneNameList) {
//            if (mappedDepthMap.containsKey(gene)) {
//                mappedDepthMap.put(gene, mappedDepthMap.get(gene).intValue() + 1);
//            } else {
//                mappedDepthMap.put(gene, 1);
//            }
//
//        }
//
//        System.out.println("bam covers " + mappedDepthMap.size() + " genes");
//
//        try {
//            BufferedWriter bw;
//            bw = YaoIOUtils.getTextGzipWriter(outFile);
//            for (Integer mappedDepth : mappedDepthMap.values()) {
//                bw.write(mappedDepth.intValue());
//                bw.newLine();
//            }
//        } catch (Exception e) {
//            e.printStackTrace();
//            System.exit(1);
//        }
//    }
// --------------------------------create position lists for UTRs,exon,gene from gff3-----------------------
    // Gene Element position collections , can be sued to next GERP element distribution Analysis
    private void callChrElementPos() {
        String gff3Dir = "/data1/home/lipeng/database/genome/wheat/annotation/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC.longesttranscript.gff3";
        String subGenome = "D";
        int geneStartPos = 0;
        int geneEndPos = 0;
        int exonStartPos = 0;
        int exonEndPos = 0;
        int threeUTRStartPos = 0;
        int threeUTREndPos = 0;
        int fiveUTRStartPos = 0;
        int fiveUTREndPos = 0;

//hash AGenome (chr1------>hash("UTR"------->pos))
        try {
            for (int i = 1; i <= 7; i++) {
                BufferedReader br;
                if (gff3Dir.endsWith("gz")) {
                    br = YaoIOUtils.getTextGzipReader(gff3Dir);
                } else {
                    br = YaoIOUtils.getTextReader(gff3Dir);
                }

                HashMap<String, ArrayList<Integer>> chrElementPosMap = new HashMap();
                ArrayList<Integer> genePosList = new ArrayList();
                ArrayList<Integer> exonPosList = new ArrayList();
                ArrayList<Integer> threeUTRPosList = new ArrayList();
                ArrayList<Integer> fiveUTRPosList = new ArrayList();

                StringBuilder chr = new StringBuilder();
                chr.append("chr").append(String.valueOf(i)).append(subGenome);
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith(chr.toString())) {
                        List<String> tList = PStringUtils.fastSplit(temp);
                        tem = tList.toArray(new String[tList.size()]);
                        if (tem[2].equals("gene")) {
                            geneStartPos = Integer.parseInt(tem[3]);
                            geneEndPos = Integer.parseInt(tem[4]);
                            for (int j = geneStartPos; j <= geneEndPos; j++) {
                                genePosList.add(j);
                            }
                        } else {
                            if (tem[2].equals("exon")) {
                                exonStartPos = Integer.parseInt(tem[3]);
                                exonEndPos = Integer.parseInt(tem[4]);
                                for (int k = exonStartPos; k <= exonEndPos; k++) {
                                    exonPosList.add(k);
                                }

                            } else {
                                if (tem[2].equals("five_prime_UTR")) {
                                    fiveUTRStartPos = Integer.parseInt(tem[3]);
                                    fiveUTREndPos = Integer.parseInt(tem[4]);
                                    for (int l = fiveUTRStartPos; l <= fiveUTREndPos; l++) {
                                        fiveUTRPosList.add(l);
                                    }
                                } else {
                                    if (tem[2].equals("three_prime_UTR")) {
                                        threeUTRStartPos = Integer.parseInt(tem[3]);
                                        threeUTREndPos = Integer.parseInt(tem[4]);
                                        for (int m = threeUTRStartPos; m <= threeUTREndPos; m++) {
                                            threeUTRPosList.add(m);
                                        }
                                    }
                                }

                            }

                        }
                    }
                }
                br.close();
                chrElementPosMap.put("gene", genePosList);
                chrElementPosMap.put("exon", exonPosList);
                chrElementPosMap.put("five_prime_UTR", fiveUTRPosList);
                chrElementPosMap.put("three_prime_UTR", threeUTRPosList);
                ElementPosMap.put(i, chrElementPosMap);
            }

        } catch (Exception e) {
            System.out.println("Error in calling " + "wheat" + subGenome + " genomic elements' Position");
            e.printStackTrace();
        }
    }

    //备选 for callChrElementPos
    private void extratLongestTranscript() {
        String gff3Dir = "/data1/home/lipeng/Hordeum_vulgare.IBSC_v2.43.gff3.gz";
        String out2Dir = "/data1/home/lipeng/Hordeum_vulgare.IBSC_v2.43.longestTrans.gff3";
        String mRNAID = null;
        int mRNALength = 0;
        ArrayList<String> infoList = new ArrayList();
        ArrayList<String> geneList = new ArrayList();
        HashMap<String, Integer> mRNAMap = new HashMap<>();
//hash AGenome (chr1------>hash("UTR"------->pos))
        try {

            BufferedReader br;
            br = YaoIOUtils.getTextGzipReader(gff3Dir);
            while ((temp = br.readLine()) != null) {

                if (!temp.startsWith("chr")) {
                    continue;
                }

                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (tem[2].contains("gene")) {
                    String[] te = tem[8].split("=");
                    if(te[1].split(";")[0].split(":")[1].startsWith("HOR")){
                    geneList.add(te[1].split(";")[0].split(":")[1]);}

                }

                infoList.add(temp);
            }
            br.close();

            String[] geneNames = geneList.toArray(new String[geneList.size()]);
            Arrays.sort(geneNames);
            System.out.println("HC has " + geneNames.length + " genes");
            String[] info = infoList.toArray(new String[infoList.size()]);
            for (int i = 0; i < info.length; i++) {
                List<String> tList = PStringUtils.fastSplit(info[i]);
                tem = tList.toArray(new String[tList.size()]);
                if (tem[2].startsWith("mRNA")) {
                    mRNALength = Integer.parseInt(tem[4]) - Integer.parseInt(tem[3]);
                    mRNAID = tem[8].split(";")[0].split(":")[1];
                    mRNAMap.put(mRNAID, mRNALength);
                }
            }

            HashMap<Integer, String> altmRNAMap = new HashMap<>();
            ArrayList<String> longestmRNAList = new ArrayList<>();
            for (String i : geneNames) {
                for (String j : mRNAMap.keySet()) {
                    if (j.startsWith(i)) {
                        altmRNAMap.put(mRNAMap.get(j), j);
                    }
                }
                ArrayList<Integer> altmRNAList = new ArrayList<>(altmRNAMap.keySet());
                Collections.sort(altmRNAList);
                
                longestmRNAList.add(altmRNAMap.get(altmRNAList.get(altmRNAList.size() - 1)));
                altmRNAMap = new HashMap<>();
            }
            System.out.println("extract " + longestmRNAList.size() + " longest transcripts");
            //write longest transcript info list
//           BufferedWriter bw;
//           
//                bw = YaoIOUtils.getTextWriter(out1Dir);
//                for(String k:longestmRNAList){
//                bw.write(k);
//                bw.newLine();
//                }
//                bw.flush();
//             bw.close();
            Collections.sort(longestmRNAList);
            BufferedWriter bw;
            bw = YaoIOUtils.getTextWriter(out2Dir);
            for (int i = 0; i < info.length; i++) {

                List<String> tList = PStringUtils.fastSplit(info[i]);
                tem = tList.toArray(new String[tList.size()]);
                if (tem[2].startsWith("gene")) {
                       if(tem[8].split("=")[1].split(";")[0].split(":")[1].startsWith("HOR")){
                    bw.write(info[i]);
                    bw.newLine();}
                } else {
                    if (tem[2].startsWith("mRNA")) {
                        mRNAID = tem[8].split(";")[0].split(":")[1];
                        if (Collections.binarySearch(longestmRNAList, mRNAID) >= 0) {
                            bw.write(info[i]);
                            bw.newLine();
                        }
                    } else {

                        mRNAID = tem[8].split(";")[0].split(":")[1];
                        if (Collections.binarySearch(longestmRNAList, mRNAID) >= 0) {
                            bw.write(info[i]);
                            bw.newLine();
                        }
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in calling " + "wheat" + subGenome + " genomic elements' Position");
            e.printStackTrace();
        }
    }

    //read in clip modified longest transcript.gff3
    public void callgeneticElementLength() {
        String gff3Dir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC.longesttranscript.gff3";
        try {
            int geneLength = 0;
            int exonLength = 0;
            int fiveUtrLength = 0;
            int threeUtrLength = 0;
            BufferedReader br;
            br = YaoIOUtils.getTextReader(gff3Dir);
            while ((temp = br.readLine()) != null) {
                if (!temp.startsWith("chr")) {
                    continue;
                }

                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (tem[0].endsWith("B")) {
                    if (tem[2].contains("three_prime_UTR")) {
                        threeUtrLength = threeUtrLength + Integer.parseInt(tem[4]) - Integer.parseInt(tem[3]);
                    }
                    if (tem[2].contains("five_prime_UTR")) {
                        fiveUtrLength = fiveUtrLength + Integer.parseInt(tem[4]) - Integer.parseInt(tem[3]);
                    }
                    if (tem[2].contains("exon")) {
                        exonLength = exonLength + Integer.parseInt(tem[4]) - Integer.parseInt(tem[3]);
                    }
                    if (tem[2].contains("gene")) {
                        geneLength = geneLength + Integer.parseInt(tem[4]) - Integer.parseInt(tem[3]);
                    }
                }
            }
            br.close();
            System.out.println("geneLength is " + geneLength + " bp");
            System.out.println("exonLength is " + exonLength + " bp");
            System.out.println("5'UTRLength is " + fiveUtrLength + " bp");
            System.out.println("3'UTRLength is " + threeUtrLength + " bp");

        } catch (Exception e) {
            System.out.println("Error in calling " + "wheat" + subGenome + " genomic elements' Position");
            e.printStackTrace();
        }

    }

//----------------------------------intragenic GERP distribution------------------------------
    //flow: count GERP>? from GERP normalized file, list GERP>? positions, binarysearch GERP>? positions in gene elenments positions list
    //output GERP numbers in each elements.
    private void intrageneGERP() {
        String subGenome = "D";
        String[] gerpFile = {"/data1/home/lipeng/result/GERP/axt/D/gerp/wheatD.chr1.gerp++",
            "/data1/home/lipeng/result/GERP/axt/D/gerp/wheatD.chr2.gerp++",
            "/data1/home/lipeng/result/GERP/axt/D/gerp/wheatD.chr3.gerp++",
            "/data1/home/lipeng/result/GERP/axt/D/gerp/wheatD.chr4.gerp++",
            "/data1/home/lipeng/result/GERP/axt/D/gerp/wheatD.chr5.gerp++",
            "/data1/home/lipeng/result/GERP/axt/D/gerp/wheatD.chr6.gerp++",
            "/data1/home/lipeng/result/GERP/axt/D/gerp/wheatD.chr7.gerp++"};
        long constraintPos0 = 0;
        long constraintPos1 = 0;
        long constraintPos2 = 0;
        HashMap<Integer, ArrayList<Integer>> gerpPosMap = new HashMap();
        ArrayList<Integer> gerpPosList = new ArrayList();
        try {
            for (int i = 0; i < gerpFile.length; i++) {
                BufferedReader br;
                if (gerpFile[i].endsWith("gz")) {
                    br = YaoIOUtils.getTextGzipReader(gerpFile[i]);
                } else {
                    br = YaoIOUtils.getTextReader(gerpFile[i]);
                }
                int pos = 1;

                while ((temp = br.readLine()) != null) {
                    List<String> tList = PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    BigDecimal gerp = new BigDecimal(tem[1]);
                    BigDecimal n = new BigDecimal("0");
                    BigDecimal m = new BigDecimal("1");
                    BigDecimal o = new BigDecimal("2");
                    if (gerp.compareTo(n) > 0) {
                        constraintPos0++;
                    }
                    if (gerp.compareTo(m) > 0) {
                        gerpPosList.add(pos);
                        constraintPos1++;
                    }
                    if (gerp.compareTo(o) > 0) {

                        constraintPos2++;
                    }
                    pos++;

                }
                gerpPosMap.put(i + 1, gerpPosList);
                br.close();
            }

            System.out.println("wheat" + subGenome + " contains " + constraintPos0 + " gerp>0 constraint Positions");
            System.out.println("wheat" + subGenome + " contains " + constraintPos1 + " gerp>1 constraint Positions");
            System.out.println("wheat" + subGenome + " contains " + constraintPos2 + " gerp>2 constraint Positions");
        } catch (Exception e) {
            System.out.println("Error in calling " + "wheat" + subGenome + " gerp Position");
            e.printStackTrace();
        }
        HashMap<Integer, HashMap<String, ArrayList<Integer>>> ElementPosMap = this.ElementPosMap;
        int constraintGenebp = 0;
        int constraintExonbp = 0;
        int constraintFiveUTRbp = 0;
        int constraintThreeUTRbp = 0;
        for (int i = 1; i <= 7; i++) {
            HashSet<Integer> genePosSet = new HashSet(ElementPosMap.get(i).get("gene"));
            HashSet<Integer> exonPosSet = new HashSet(ElementPosMap.get(i).get("exon"));
            HashSet<Integer> threeUTRPosSet = new HashSet(ElementPosMap.get(i).get("three_prime_UTR"));
            HashSet<Integer> fiveUTRPosSet = new HashSet(ElementPosMap.get(i).get("five_prime_UTR"));
            ArrayList<Integer> genePosList = new ArrayList();
            ArrayList<Integer> exonPosList = new ArrayList();
            ArrayList<Integer> threeUTRPosList = new ArrayList();
            ArrayList<Integer> fiveUTRPosList = new ArrayList();
            genePosList.addAll(genePosSet);
            exonPosList.addAll(exonPosSet);
            fiveUTRPosList.addAll(fiveUTRPosSet);
            threeUTRPosList.addAll(threeUTRPosSet);
            Collections.sort(genePosList);
            Collections.sort(exonPosList);
            Collections.sort(fiveUTRPosList);
            Collections.sort(threeUTRPosList);
            for (int j : gerpPosMap.get(i)) {
                if (Collections.binarySearch(genePosList, j) >= 0) {
                    constraintGenebp++;
                }
                if (Collections.binarySearch(exonPosList, j) >= 0) {
                    constraintExonbp++;
                }
                if (Collections.binarySearch(fiveUTRPosList, j) >= 0) {
                    constraintFiveUTRbp++;
                }
                if (Collections.binarySearch(threeUTRPosList, j) >= 0) {
                    constraintThreeUTRbp++;
                }

            }
        }
        System.out.println("gerp>1 covers " + constraintGenebp + " bp of gene");
        System.out.println("gerp>1 covers " + constraintExonbp + " bp of exon");
        System.out.println("gerp>1 covers " + constraintFiveUTRbp + " bp of fiveUTR");
        System.out.println("gerp>1 covers " + constraintThreeUTRbp + " bp of threeUTR");
    }

//------------------------------------count and write GERP>? from GERP.rates raw data file-----------------------
    public void callGERPConstraintNum(String subGenome) {
        String GERPDir = "/data1/home/lipeng/result/GERP/axt/D/gerp/Triticeae.mfa.jgerp";
        try {
            BufferedReader br;

            if (GERPDir.endsWith("gz")) {
                br = YaoIOUtils.getTextGzipReader(GERPDir);

            } else {
                br = YaoIOUtils.getTextReader(GERPDir);

            }
            long constraintNum = 0;
            while ((temp = br.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                BigDecimal gerp = new BigDecimal(tem[1]);
                BigDecimal n = new BigDecimal("0");
                if (gerp.compareTo(n) > 0) {
                    constraintNum++;
                }

            }
            System.out.println("wheat " + subGenome + "contains " + constraintNum + " GERP>0");
            br.close();
        } catch (Exception e) {
            System.out.println("Error in calling sample GERP!");
            e.printStackTrace();
        }

    }

    public void parseGerpelem() {
        String GERPDir = "/data1/home/lipeng/result/GERP/axt/A/gerp/wheatA.chr1.gerp++.elems";
        try {
            BufferedReader br;

            if (GERPDir.endsWith("gz")) {
                br = YaoIOUtils.getTextGzipReader(GERPDir);

            } else {
                br = YaoIOUtils.getTextReader(GERPDir);

            }

            ArrayList<Integer> gerpelemPos = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                int startPos = Integer.parseInt(tem[1]);
                int endPos = Integer.parseInt(tem[2]);
                for (int i = startPos; i <= endPos; i++) {
                    gerpelemPos.add(i);
                }
            }
            ArrayList<Integer> genePosList = this.genePosList;

            long pos = 0;
            for (int j : gerpelemPos) {
                if (Collections.binarySearch(genePosList, j) >= 0) {
                    pos++;
                }
            }
            System.out.println("gerpelem--wheat chr 1 constarint region contains " + gerpelemPos.size() + " bp");
            System.out.println("gerpelem--wheat chr 1 constarint region contains " + pos + " bp genes' region");
            br.close();
        } catch (Exception e) {
            System.out.println("Error in calling sample GERP!");
            e.printStackTrace();
        }

    }

//    private void callMappedDepth(String depthFile, String subGenome) {
//        HashMap<Integer, HashMap<String, ArrayList<Integer>>> ElementPosMap = this.ElementPosMap;
//        //call UTR,gene,intron,exon,reads mapped depth
//
//        BufferedReader br;
//        if (depthFile.endsWith("gz")) {
//            br = YaoIOUtils.getTextGzipReader(depthFile);
//        } else {
//            br = YaoIOUtils.getTextReader(depthFile);
//        }
//        int pos = 0;
//        long geneDepth = 0;
//        long fiveUTRDepth = 0;
//        long threeUTRDepth = 0;
//        long exonDepth = 0;
//        try {
//            while ((temp = br.readLine()) != null) {
//
//                for (int i = 1; i <= 7; i++) {
//                    HashSet<Integer> genePosSet = new HashSet(ElementPosMap.get(i).get("gene"));
//                    HashSet<Integer> exonPosSet = new HashSet(ElementPosMap.get(i).get("exon"));
//                    HashSet<Integer> threeUTRPosSet = new HashSet(ElementPosMap.get(i).get("three_prime_UTR"));
//                    HashSet<Integer> fiveUTRPosSet = new HashSet(ElementPosMap.get(i).get("five_prime_UTR"));
//                    ArrayList<Integer> genePosList = new ArrayList();
//                    ArrayList<Integer> exonPosList = new ArrayList();
//                    ArrayList<Integer> threeUTRPosList = new ArrayList();
//                    ArrayList<Integer> fiveUTRPosList = new ArrayList();
//                    genePosList.addAll(genePosSet);
//                    exonPosList.addAll(exonPosSet);
//                    fiveUTRPosList.addAll(fiveUTRPosSet);
//                    threeUTRPosList.addAll(threeUTRPosSet);
//                    Collections.sort(genePosList);
//                    Collections.sort(exonPosList);
//                    Collections.sort(fiveUTRPosList);
//                    Collections.sort(threeUTRPosList);
//                    StringBuilder chr = new StringBuilder();
//                    chr.append("chr").append(String.valueOf(i)).append(subGenome);
//                    if (temp.startsWith(chr.toString())) {
//                        List<String> tList = PStringUtils.fastSplit(temp);
//                        tem = tList.toArray(new String[tList.size()]);
//                        pos = Integer.parseInt(tem[1]);
//                        if (Collections.binarySearch(genePosList, pos) >= 0) {
//                            geneDepth = geneDepth + Integer.parseInt(tem[2]);
//                        }
//                        if (Collections.binarySearch(exonPosList, pos) >= 0) {
//                            exonDepth = exonDepth + Integer.parseInt(tem[2]);
//                        }
//                        if (Collections.binarySearch(fiveUTRPosList, pos) >= 0) {
//                            fiveUTRDepth = fiveUTRDepth + Integer.parseInt(tem[2]);
//                        }
//                        if (Collections.binarySearch(threeUTRPosList, pos) >= 0) {
//                            threeUTRDepth = threeUTRDepth + Integer.parseInt(tem[2]);
//                        }
//                    }
//                }
//            }
//            br.close();
//            System.out.println("--------------wheat" + subGenome + " genome element mapped bases depth counting--------------");
//            System.out.println("gene mapped bases number:" + geneDepth);
//            System.out.println("exon mapped bases number:" + exonDepth);
//            System.out.println("fiveUTR mapped bases number:" + fiveUTRDepth);
//            System.out.println("threeUTR mapped bases number:" + threeUTRDepth);
//        } catch (Exception e) {
//            System.out.println("Error in calling " + "wheat" + subGenome + " genome element mapped bases depth");
//            e.printStackTrace();
//        }
//    }
    //---------------------------------window size distribution of GERP------------------------------------  
    // used for follwing whole genome GERP window distribution line chart; workflow:add all ?kb window GERP to sum score and write these sum socres in a column. 
    public void GERPWindowDistribution() {
        String GERPDir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/Qgene20k.gerp";
        String outFile = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/Qgene.1kbwindow";
        try {
            BufferedReader br;
            BufferedWriter bw;
            if (GERPDir.endsWith("gz")) {
                br = YaoIOUtils.getTextGzipReader(GERPDir);
                bw = YaoIOUtils.getTextGzipWriter(outFile);
            } else {
                br = YaoIOUtils.getTextReader(GERPDir);
                bw = YaoIOUtils.getTextWriter(outFile);
            }
            int windowSize = 1;
            ArrayList<String> GERPsumList = new ArrayList();

            float GERPsum = 0;
            DecimalFormat df = new DecimalFormat("0.000");
            while ((temp = br.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                GERPsum = GERPsum + Float.parseFloat(tem[1].toString());

                if (windowSize == 1000) {
                    String stringGERPsum = df.format(GERPsum);
                    GERPsumList.add(stringGERPsum);
                    windowSize = 0;
                    GERPsum = 0;
                }
                windowSize++;
            }
//provide GERP lines is not x times of windowSize,add another list.add to get last one GERP sum
            String stringGERPsum = df.format(GERPsum);
            GERPsumList.add(stringGERPsum);
            for (String i : GERPsumList) {
                bw.write(i.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            System.out.println("Error in calling " + "wheat" + subGenome + " genome element mapped bases depth");
            e.printStackTrace();
        }
    }

    public void calculateMAF() {

    }

    //plot sift,gerp, synonymous and nonsynonymous
    public void sFSofDeletrious() {
        String transDir="/data1/home/lipeng/database/genome/wheat/annotation/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC.longesttranscript.gff3";
        String siftDir = "/data1/home/lipeng/result/sift/chr1A_SIFTannotations.xls";
       // String siftDir="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/chr1A_SIFTannotations.xls";
        String GERPDir = "/data1/home/lipeng/result/GERP/axt/A/gerp/wheatA.chr1.gerp++";
        //String GERPDir="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/gerp.txt";
        String afDir = "/data1/home/lipeng/result/sift/chr1A.MAF";
      //  String afDir="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.AF";
        String outFile = "/data1/home/lipeng/result/sift/longSFS.txt";
       // String outFile="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/out.txt";
        try {
             BufferedReader br0;
              br0 = YaoIOUtils.getTextReader(transDir);
              ArrayList<String> transList=new ArrayList();
              br0.readLine();
                while ((temp = br0.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if(tem[2].startsWith("mRNA")){
                transList.add(tem[8].split(";")[0].substring(3)); 
                }
 
                }
                Collections.sort(transList);
                br0.close();
            BufferedReader br1;
            BufferedWriter bw;

            br1 = YaoIOUtils.getTextReader(siftDir);
            bw = YaoIOUtils.getTextWriter(outFile);
            Set<Integer> siftPos = new HashSet();
            br1.readLine();
            String varType1 = "SYNONYMOUS";
            float siftScore2 = 0;
            String varType2 = null;
            float siftScore1 = 2.61f;
            int Pos1 = 1160416;
            int Pos2 = 0;
            HashMap<Integer, Float> siftScore = new HashMap<>();
            HashMap<Integer, String> varType = new HashMap();
            while ((temp = br1.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (tem[13].startsWith("NA"))  continue;
                
            if(Collections.binarySearch(transList,tem[4])>=0){
            siftScore.put(Integer.parseInt(tem[1]), Float.parseFloat(tem[12].toString()));
            varType.put(Integer.parseInt(tem[1]), tem[8]);
             siftPos.add(Integer.parseInt(tem[1]));
            
            }
            
            }
//                Pos2 = Integer.parseInt(tem[1]);
//                siftScore2 = Float.parseFloat(tem[12].toString());
//                varType2 = tem[8];
//                siftPos.add(Integer.parseInt(tem[1]));
//                if (Pos2 == Pos1) {
//                    if (varType2.startsWith("NONSYNONYMOUS") || varType1.startsWith("NONSYNONYMOUS")) {
//                        varType1 = "NONSYNONYMOUS";
//                    }
//                    if (siftScore2 <= siftScore1) {
//                        siftScore1 = siftScore2;
//                    }
//                } else {
//                    siftScore.put(Pos1, siftScore1);
//                    varType.put(Pos1, varType1);
//                    Pos1 = Pos2;
//                    varType1 = varType2;
//                    siftScore1 = siftScore2;
//                }
//            }
//            siftScore.put(Pos2, siftScore2);
//            varType.put(Pos2, varType2);

            br1.close();
            System.out.println("sift parsing : ok");
            ArrayList<Integer> siftPosList = new ArrayList<>(siftPos);
            Collections.sort(siftPosList);
            
            BufferedReader br2;
            br2 = YaoIOUtils.getTextReader(GERPDir);
            float GERP = 0;
            int Pos = 0;
            int Num = 0;
            HashMap<Integer, Float> gerpmap = new HashMap();
            while ((temp = br2.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                GERP = Float.parseFloat(tem[1].toString());
                Pos++;
                if (Pos == siftPosList.get(Num)) {
                    gerpmap.put(Pos, GERP);
                    Num++;
                    if (Num >= siftPosList.size()) {
                        break;
                    }
                }
            }
            br2.close();
            System.out.println("GERP parsing : ok");
            HashMap<Integer, String> afmap = new HashMap();

            BufferedReader br3;
            br3 = YaoIOUtils.getTextReader(afDir);
            br3.readLine();
            br3.readLine();
            ArrayList<Integer> afList = new ArrayList();
            while ((temp = br3.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                afList.add(Integer.parseInt(tem[0]));
                afmap.put(Integer.parseInt(tem[0]), tem[1]);
            }
            Collections.sort(afList);

            br3.close();
            System.out.println("AF   parsing: ok");
            // SNPs in CDS (with sift score) will be parse
            bw.write("Pos");
            bw.write("\t");
            bw.write("GERP");
            bw.write("\t");
            bw.write("variant_Type");
            bw.write("\t");
            bw.write("SIFT");
            bw.write("\t");
            bw.write("MAF");
            bw.newLine();

            for (int j : siftPosList) {
                if (Collections.binarySearch(afList, j) >= 0) {
                    bw.write(String.valueOf(j));
                    bw.write("\t");
                    bw.write(String.valueOf(gerpmap.get(j)));
                    bw.write("\t");
                    bw.write(varType.get(j));
                    bw.write("\t");
                    bw.write(String.valueOf(siftScore.get(j)));
                    bw.write("\t");
                    bw.write(String.valueOf(afmap.get(j)));
                    bw.newLine();

                }
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            System.out.println("Error in calling " + "wheat" + subGenome + " genome element mapped bases depth");
            e.printStackTrace();
        }
    }

    public void deleteriousstat() {
        String infile = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/longSFS.txt";
        try {
            BufferedReader br;
            br = YaoIOUtils.getTextReader(infile);
            br.readLine();
            int synNum = 0;
            int nsynNum = 0;
            float synSum = 0;
            float nsynSum = 0;
            float sumGERP = 0;
            float jzMAF = 0;
            double[] MAF = {0, 0.02, 0.04, 0.06, 0.08, 0.10,
                0.12, 0.14, 0.16, 0.18, 0.20,
                0.22, 0.24, 0.26, 0.28, 0.30,
                0.32, 0.34, 0.36, 0.38, 0.40,
                0.42, 0.44, 0.46, 0.48, 0.50,};
            ArrayList<ArrayList<Float>> GERPList = new ArrayList<>();
            ArrayList <String> lines=new ArrayList<>();
            int[]nsynGERP=new int[25];
            for (int l = 0; l <= 24; l++) {
            nsynGERP[l]=0;
            }
             int[]synGERP=new int[25];
            for (int l = 0; l <= 24; l++) {
           synGERP[l]=0;
            }
            for (int l = 0; l <= 24; l++) {
                ArrayList<Float> emptyGERP = new ArrayList<>();
                GERPList.add(emptyGERP);
            }
            while ((temp = br.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                for (int i = 0; i < 25; i++) {
//                    if (Float.parseFloat(tem[4]) > 0.5) {
//                        jzMAF = 1 - Float.parseFloat(tem[4]);
//                    } else {
//                        jzMAF = Float.parseFloat(tem[4]);
//                    }
jzMAF = Float.parseFloat(tem[4]);
                    if (jzMAF <= MAF[i + 1] && jzMAF > MAF[i]) {
                        GERPList.get(i).add(Float.parseFloat(tem[1]));
                    if(jzMAF<=0.34&&jzMAF>0.32){lines.add(temp);}}}}
//                      if(tem[2].startsWith("SYNONYMOUS")){
//                          synGERP[i]++;
//                      }
//                       if(tem[2].startsWith("NONSYNONYMOUS")){
//                     if (Float.parseFloat(tem[1])>0 && Float.parseFloat(tem[3])<=0.05){
//                     nsynGERP[i]++;
//                     }}}}}
            
            br.close();
//            for (int j = 0; j <= 24; j++) {
//             System.out.println(nsynGERP[j]+"\t"+synGERP[j]);
//            }
           //######
            for (int j = 0; j <= 24; j++) {
//                for (float k : GERPList.get(j)) {
////                    sumGERP = sumGERP + k;
//
////                }
//                System.out.println(sumGERP);
//                sumGERP = 0;
System.out.println(GERPList.get(j).size());
            }
            for (String j:lines){System.out.println(j);}
//#### SYN and NSYN GERP distribution
//       if(Float.parseFloat(tem[4])>0&&Float.parseFloat(tem[4])<=0.05){
//       if(tem[2].startsWith("SYNONYMOUS")){SYNONYMOUS++;}
//       }
//       if(tem[2].startsWith("SYNONYMOUS")){synNum++;synSum=synSum+Float.parseFloat(tem[1]);}else{
//       nsynNum++;nsynSum=nsynSum+Float.parseFloat(tem[1]);
//       }
//       }
//       System.out.println("SYNONYMOUS: "+synSum/synNum);
//       System.out.println("NONSYNONYMOUS: "+nsynSum/nsynNum);
        } catch (Exception e) {
            System.out.println("Error in calling " + "wheat" + subGenome + " genome element mapped bases depth");
            e.printStackTrace();
        }
    }
    
    public void figMAF(){
     String infile = "/data1/home/lipeng/result/sift/chr1A.filt.MafMiss.vcf";
     String outfile = "/data1/home/lipeng/result/sift/chr1A.MafMiss.MAF";
     try{
          BufferedReader br = YaoIOUtils.getTextReader(infile);
          BufferedWriter bw = YaoIOUtils.getTextWriter(outfile);
          
          bw.write("CHR1A-POS");bw.write("\t");bw.write("MAF");bw.newLine();bw.write("author:clip lee");bw.newLine();
     while ((temp = br.readLine()) != null) {
     if(!temp.startsWith("1A"))continue;
     int DP=0; 
      int AD1=0;
     int AD2=0;
     int AD3=0;
     List<String> tList = PStringUtils.fastSplit(temp);
     tem = tList.toArray(new String[tList.size()]);
     if(tem[4].length()>1){
        bw.write(tem[1]);bw.write("\t");
     for(int i=9;i<tem.length;i++){
         if(tem[i].startsWith("./."))continue;
  String te[]= tem[i].split(":")[1].split(",");
  for(int j=0;j<te.length;j++){
      DP=DP+Integer.parseInt(te[j]);
  }
  AD1=AD1+Integer.parseInt(te[0]);
  AD2=AD2+Integer.parseInt(te[1]);
  AD3=AD3+Integer.parseInt(te[2]);
     }
  ArrayList<Integer> ADList =new ArrayList();
  ADList.add(AD1);  ADList.add(AD2);  ADList.add(AD3);
  Collections.sort(ADList);
 
      bw.write(String.valueOf((double) Math.round(1000*ADList.get(1)/DP)/1000));
      bw.newLine();

     }else{
       bw.write(tem[1]);bw.write("\t");
     for(int i=9;i<tem.length;i++){
         if(tem[i].startsWith("./."))continue;
  String te[]= tem[i].split(":")[1].split(",");
  for(int j=0;j<te.length;j++){
      DP=DP+Integer.parseInt(te[j]);
  }
  AD1=AD1+Integer.parseInt(te[0]);
  AD2=AD2+Integer.parseInt(te[1]);}
  
  if(AD1>=AD2){
      bw.write(String.valueOf((double) Math.round(1000*AD2/DP)/1000));
      bw.newLine();
  }else{
      bw.write(String.valueOf((double) Math.round(1000*AD1/DP)/1000));
      bw.newLine();}
     }}
     br.close();
    bw.flush();
    bw.close();
      } catch (Exception e) {
            System.out.println("Error in calling " + "wheat" + subGenome + " genome element mapped bases depth");
            e.printStackTrace();
        }
     }
    
    

// ----------------------------------parse parameters--------------------------------
    //this function is used to get parameters, usually use it but sometimes not
    private void parseParameters(String infileS) {
        List<String> pLineList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            boolean ifOut = false;
            if (!(temp = br.readLine()).equals("calling maf cover gene")) {
                ifOut = true;
            }
            if (!(temp = br.readLine()).equals("Author: clip")) {
                ifOut = true;
            }
            if (!(temp = br.readLine()).equals("Email: lipenkang@163.com")) {
                ifOut = true;
            }
            if (!(temp = br.readLine()).equals("Homepage: http://plantgeneticslab.weebly.com/")) {
                ifOut = true;
            }
            if (ifOut) {
                System.out.println("Thanks for using TEP.");
                System.out.println("Please keep the authorship in the parameter file. Program stops.");
                System.exit(0);
            }
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    continue;
                }
                if (temp.isEmpty()) {
                    continue;
                }
                pLineList.add(temp);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.gff3Dir = pLineList.get(0);
        this.depthFile = pLineList.get(2);
        this.subGenome = pLineList.get(1);
        this.outFile = pLineList.get(4);
        this.samFile = pLineList.get(3);

    }

    public static void main(String[] args) {
        new wheatGeneQual();

    }

}
