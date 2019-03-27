/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

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
       // this.callGeneLength();//this function input is standard IWSC gff3 ,not gff3 modified by fei
        this.indexGene();     //index genes
        this.parseGerpelem();
      // this.callMappedGenePos();
        // this.callGeneMappedDepth(outFile);           //need update .......
        // this.callChrElementPos(gff3Dir, subGenome);  //create position lists for UTRs,exon,gene from gff3
        // this.callMappedDepth(depthFile, subGenome); // need update.......
        // this.intrageneGERP(subGenome);          // (require callChrElementPos) GERP distribution in intragenetic regions
        //this.callGERPConstraintNum(subGenome);    //call GERP >? regions
        //this.GERPWindowDistribution();            //GERP ?bp window distribution of chromosomes
    }
    String gff3Dir = null;
    String subGenome = null;
    String outFile = null;
    String samFile = null;
    String depthFile = null;
    String temp = null;
    String[] tem = null;
    TIntArrayList mappedPos = new TIntArrayList();
    HashMap<Integer, ArrayList<Integer>> mappedReadsPosMap = new HashMap<>();
    HashMap<Integer, String> chrGeneIndexMap = new HashMap<>();
    HashMap<Integer, HashMap<Integer, String>> allGeneIndexMap = new HashMap<>();
    ArrayList<Integer> genePosList = new ArrayList();
    HashMap<Integer, ArrayList<Integer>> allGenePosMap = new HashMap<>();
    HashMap<String, Integer> mappedDepthMap = new HashMap<>();
    HashMap<Integer, HashMap<String, ArrayList<Integer>>> ElementPosMap = new HashMap();
    ArrayList<String> readsNameList = new ArrayList<>();

//--------------------------------------call gene or intrageneic elements length or numbers from gff3------------------------
    //used to following hist plot gene length
    private void callGeneLength() {
      // String gff3Dir = "/data1/home/lipeng/database/genome/urartu/WheatTu.annotation/WheatTu.gene.gff";
        String gff3Dir ="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.gff3";
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

                    if (tem[2].equals("gene") && tem[0].length()==3) {
                        StartPos = Integer.parseInt(tem[3]);
                        EndPos = Integer.parseInt(tem[4]);
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
//--------------------------indexGene positions--------------------------
    //create index of all gene Pos in GFF3

    private  ArrayList<Integer> indexGene() {
      // String gff3Dir = "/data1/home/lipeng/database/genome/urartu/WheatTu.annotation/WheatTu.gene.gff";
        //      String gff3Dir ="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.gff3";
 String gff3Dir="/data1/home/lipeng/database/genome/wheat/annotation/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gff3";
       BufferedReader br;
        if (gff3Dir.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(gff3Dir);
        } else {
            br = YaoIOUtils.getTextReader(gff3Dir);
        }
        int geneStartPos = 0;
        int geneEndPos = 0;
        String geneInfo = null;

       // ArrayList<Integer>[] genePosListArr = new ArrayList[7];
//        for (int i = 0; i <= 6; i++) {
//            genePosListArr[i] = new ArrayList<Integer>();
//        }
        // HashMap[] geneIndexArr = new HashMap[7];
        
        try {

            while ((temp = br.readLine()) != null) {
//                if (temp.startsWith("TuUngrouped")) {
//                    continue;
//                }
 if (temp.startsWith("chr1A")) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
//                    StringBuilder chr = new StringBuilder();
//                    chr.append(String.valueOf(i)).append("A");
                //int chr = Integer.valueOf(String.valueOf(tem[0].charAt(3)));

                if (tem[2].equals("gene")) {
                    if(geneStartPos>=1000){
                    geneStartPos = Integer.parseInt(tem[3])-10000;
                    geneEndPos = Integer.parseInt(tem[4])+10000;}else{
                    geneStartPos = 0;
                     geneEndPos = Integer.parseInt(tem[4])+10000;
                    }
                   // geneInfo = KStringUtils.fastSplitSemicolon(tem[8]).get(0).substring(3);
                    // hash all gene Position in a chr ------>gene name
                    for (int j = geneStartPos; j <= geneEndPos; j++) {
                        //  geneIndexArr[chr].put(j, geneInfo);
                       // genePosListArr[chr - 1].add(j);
                       genePosList.add(j);
                    }
                }
            }}
            
            br.close();
//            for (int i = 1; i < 8; i++) {
//                allGenePosMap.put(i, genePosListArr[i - 1]);
//                //  allGeneIndexMap.put(i, geneIndexArr[i]);
//            }
//                genePosList = new ArrayList<Integer>();
//                chrGeneIndexMap = new HashMap<Integer, String>();
System.out.println("indexing genes in urartu genome. ok!");
        } catch (Exception e) {
            System.out.println("Error in indexing genes in urartu genome");
            e.printStackTrace();
        }
        return genePosList;
    }
//-----------------------------parse minimap2.PAF ------------------------------------------
    //collect mapped contigs pos 

    public void callMappedGenePos() {
       //String pafFile = "/data1/home/lipeng/result/PacBioAnalysis/urartuQC1/readschain.kbmap";
         String pafFile ="/Users/kanglipeng/Desktop/reads.kbmap";
        BufferedReader br;
        if (pafFile.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(pafFile);
        } else {
            br = YaoIOUtils.getTextReader(pafFile);
        }
        int mappedBlockStart = 0;
        int mappedBlockEnd = 0;
        long contigsSumLength = 0;
//        ArrayList[] mappedPosArr = new ArrayList[7];
//        for (int i = 0; i <= 6; i++) {
//            mappedPosArr[i] = new ArrayList<Integer>();
//        }
ArrayList<Integer> mappedPosList=new ArrayList<>();
        try {
            while ((temp = br.readLine()) != null) {

                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (tem[5].equals("Tu1")) {
                    
                mappedBlockStart = Integer.parseInt(tem[8]);
                mappedBlockEnd = Integer.parseInt(tem[9]);
                contigsSumLength = contigsSumLength + mappedBlockEnd - mappedBlockStart;
              //  int chr = Integer.valueOf(String.valueOf(tem[5].charAt(2)));
                for (int j = mappedBlockStart; j <= mappedBlockEnd; j++) {
                    //mappedPosArr[chr - 1].add(j);
                    mappedPosList.add(j);
                }

            }}
            System.out.println("indexing wtdbg2-contigs. ok!");
            br.close();
            //create mapped contigs pos lists for each chr 
//            for (int i = 1; i <= 7; i++) {
//                Collections.sort(mappedPosArr[i - 1]);
//                mappedReadsPosMap.put(i, mappedPosArr[i - 1]);
//            }

          System.out.println(" wtdbg2-contigs index map constructing. ok!");

            //start find intersection between contigs and reference

            HashMap<Integer, ArrayList<Integer>> allGenePosMap = this.allGenePosMap;
            for (int j = 1; j <= 7; j++) {
               ArrayList <Integer>interSectionPos = new ArrayList<>();
               interSectionPos= range.getCommonElements(mappedPosList,allGenePosMap.get(j));
               
//                for (int k : allGenePosMap.get(j)) {
//                    if (Collections.binarySearch(mappedReadsPosMap.get(j), k) >= 0) {
//                        interSectionPos.add(k);
//                    }
//                }
                System.out.println("Urartu chr " + j + " was covered " + interSectionPos.size() + " bp by wtdbg2-contigs");
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
    private void callChrElementPos(String gff3Dir, String subGenome) {
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
//----------------------------------intragenic GERP distribution------------------------------
    //flow: count GERP>? from GERP normalized file, list GERP>? positions, binarysearch GERP>? positions in gene elenments positions list
    //output GERP numbers in each elements.

    private void intrageneGERP(String subGenome) {
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
    public void parseGerpelem(){
     String GERPDir = "/data1/home/lipeng/result/GERP/axt/A/gerp/wheatA.chr1.gerp++.elems";
        try {
            BufferedReader br;

            if (GERPDir.endsWith("gz")) {
                br = YaoIOUtils.getTextGzipReader(GERPDir);

            } else {
                br = YaoIOUtils.getTextReader(GERPDir);

            }
            
            ArrayList<Integer>gerpelemPos=new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                int startPos=Integer.parseInt(tem[1]);
                 int endPos=Integer.parseInt(tem[2]);
                 for(int i=startPos;i<=endPos;i++){
                 gerpelemPos.add(i);
                 }
            }
              ArrayList<Integer> genePosList= this.genePosList;
             
              long pos=0;
            for(int j: gerpelemPos){
           if(Collections.binarySearch(genePosList,j)>=0){
            pos++;
            }}
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
