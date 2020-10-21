/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import com.koloboke.collect.map.hash.HashIntFloatMap;
import com.koloboke.collect.map.hash.HashIntFloatMaps;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import utils.IOUtils;
import utils.PStringUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author kanglipeng
 */
//this package contain tools for format convert
public class FormatConvertTools {

    String temp = null;
    String[] tem = null;
    Options options = new Options();
    String faiDir = null;
    String mafDir = null;
    String outDir = null;
    String specie = null;
    String subGenome = null;

    public FormatConvertTools(String[] args) {
        this.createOptions();
        this.retrieveParameters(args);
        //this.simpleMsa();
        // this.est(mafDir, outDir,subGenome);
        //this.dfe(mafDir, outDir,subGenome);
        //  this.focal(mafDir, outDir,subGenome);
        //   this.mafToFa();
     //   this.mafTophy();
        //this.gerpReformat();
        this.mafRevise(mafDir,outDir,specie,faiDir);
        //  this.tata(mafDir, outDir);
    }
//---------------------msa to non-deletion msa--------------------------------------------------
    //convert gtf from IWGSC to sift4g input gtf
    //sift-gtf requires 'gene_biotype "protein_coding"' in ninth column of normal gtf

    public void simpleMsa() {
        BufferedReader br = IOUtils.getTextReader(mafDir);
        BufferedWriter bw = YaoIOUtils.getTextWriter(outDir);
        ArrayList<Integer> deletionPos = new ArrayList<>();
        String seq = null;
        try {
            while ((temp = br.readLine()) != null) {
                if (!temp.startsWith("a") && !temp.startsWith("s")) {
                    bw.write(temp);
                    bw.newLine();
                    continue;
                }
                if (temp.startsWith("a")) {
                    bw.write(temp);
                    bw.newLine();
                    deletionPos = new ArrayList<>();
                    continue;
                }
                tem = KStringUtils.mafSplit(temp);
                if (temp.startsWith("s tra")) {

                    for (int i = 0; i < tem[6].length(); i++) {
                        if (tem[6].charAt(i) == '-') {
                            deletionPos.add(i);
                        }
                    }
                    for (int j = 0; j < 6; j++) {
                        bw.write(tem[j]);
                        bw.write("\t");
                    }
                    seq = tem[6].replace("-", "");
                    bw.write(seq);
                    bw.newLine();
                } else {

                    seq = tem[6];
                    int rmInsertionBps = 0;
                    StringBuilder temSeq = new StringBuilder();
                    for (int k = 0; k < seq.length(); k++) {
                        if (deletionPos.contains(k)) {
                            if (seq.charAt(k) == '-') {
                                rmInsertionBps++;
                            }
                            continue;
                        }
                        temSeq.append(seq.charAt(k));
                    }
                    for (int j = 0; j < 3; j++) {
                        bw.write(tem[j]);
                        bw.write("\t");
                    }
                    bw.write(String.valueOf(Integer.parseInt(tem[3]) - deletionPos.size() + rmInsertionBps));
                    bw.write("\t");
                    for (int j = 4; j < 6; j++) {
                        bw.write(tem[j]);
                        bw.write("\t");
                    }
                    bw.write(temSeq.toString());
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in simplifying msa! please try again! Come on !");
            e.printStackTrace();
        }
        System.out.println("Congratulations! remove insertion and deletion in Ref of msa sucessfully!");
    }

    //mafToFa v1.1 : convert maf format to fa format. pay attention! this convertor only used to bring up fa format for GERP++
    public void mafToFa() {
        String maf = "/data2/lipeng/msa/multiz/a.4d.maf";
        String out = "/data2/lipeng/msa/multiz/a.4d.fa";
        BufferedReader br = YaoIOUtils.getTextReader(maf);
        BufferedWriter bw = YaoIOUtils.getTextWriter(out);
         if (maf.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(maf);
        } else {
            br = YaoIOUtils.getTextReader(maf);
        }

        String[] species = {"traes", "aetau", "trura", "hospo", "hovul", "brdis", "phedu", "leper", "orbra", "orpun", "ormer", "orlon", "orglu", "orbar", "orgla", "orniv", "orind", "orruf", "orjap", "ertef", "ortho", "sobic", "zemay", "seita", "pahal", "pavir"};
        HashMap<String, String> blockInfo = new HashMap();
        ArrayList<String>[] specieSeq = new ArrayList[species.length];
        for (int l = 0; l < species.length; l++) {
            specieSeq[l] = new ArrayList<String>();
        }
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("s")) {
                    tem = KStringUtils.mafSplit(temp);
                    blockInfo.put(tem[1].split("\\.")[0], tem[6]);
                }

            }
            br.close();
            bw.write("26 ");
            bw.write(out);
            bw.newLine();
            for (int k = 0; k < specieSeq.length; k++) {

                bw.write(">" + species[k]);
                bw.newLine();
                for (String m : specieSeq[k]) {
                    bw.write(m);
                    bw.newLine();
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in convert filtered maf to phy! please try again! Come on !");
            e.printStackTrace();
        }
        System.out.println("Congratulations! filtered maf converts to phy sucessfully!");
    }

    //usage java -jar maf2phy.jar -m maf -o out -s speciesA,speciesB,..speciesN(maf_ordered )
    public void mafTophy() {
        BufferedReader br = YaoIOUtils.getTextReader(mafDir);
        BufferedWriter bw = YaoIOUtils.getTextWriter(outDir);
        String species[] = KStringUtils.fastSplitComma(specie);
        int specieNum = 0;
        String fourDbase = null;
        // String[] species = {"traes", "aetau", "trura", "hospo", "hovul", "brdis", "phedu", "leper", "orbra", "orpun", "ormer", "orlon", "orglu", "orbar", "orgla", "orniv", "orind", "orruf", "orjap", "ertef", "ercur", "ortho", "sobic", "zemay", "seita", "eccru", "pahal", "pavir", "ancom"};
        // String[] species = {"traes","secer", "hospo", "hovul", "brdis", "phedu", "leper", "orbra", "orpun", "ormer", "orlon", "orglu", "orbar", "orgla", "orniv", "orind", "orruf", "orjap", "ertef", "ercur", "ortho", "sobic","coaqu", "zemay", "seita", "eccru", "pahal", "pavir", "ancom"};
        ArrayList<String>[] specieSeq = new ArrayList[species.length];
        for (int l = 0; l < species.length; l++) {
            specieSeq[l] = new ArrayList<String>();
        }
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("a")) {
                    specieNum = 0;
                }
                if (temp.startsWith("s")) {
                    tem = KStringUtils.mafSplit(temp);
                    //mafFilter & mafsInRegion retain deletion and insertion "-" appended column
                    if (tem[6].length() > 1) {
                        fourDbase = tem[6].substring(0, 1);
                    } else {
                        fourDbase = tem[6];
                    }
                    specieSeq[specieNum].add(fourDbase);
                    specieNum++;
                }
            }
            br.close();
            bw.write(String.valueOf(species.length) + " ");
            bw.write(String.valueOf(specieSeq[0].size()));
            bw.newLine();
            for (int k = 0; k < specieSeq.length; k++) {
                bw.write(species[k]);
                bw.write("\t");
                for (String m : specieSeq[k]) {
                    bw.write(m);
                }
                bw.newLine();
                System.out.println(species[k] + "\t" + String.valueOf(specieSeq[0].size()));
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            System.out.println("Error in convert maf to phy! please try again! Come on !");
            e.printStackTrace();
        }
        System.out.println("Congratulations! MAF converts to phy sucessfully!");
    }
    
    
    

    //gerpNorm: convert GERP.rates to phylop standard output format
    //example:
// GERP.rates                       wig format(1-based)
//0       0       ---->          fixedStep chrom=chr1A start=1 step=1
//0       0                      0 
//0       0       ---->          0
//0.678   0.678                  0
//                ---->          0.678
    //-f x.bed -m x.rates -o x.out -g A 
    public void gerpReformat() {

        BufferedReader brg = YaoIOUtils.getTextReader(mafDir);
        BufferedReader brb = YaoIOUtils.getTextReader(faiDir);//mafRanges bed
        BufferedWriter bw = YaoIOUtils.getTextWriter(outDir);
        try {
            // HashIntFloatMap gerpScoreMap = HashIntFloatMaps.newMutableMap();
            float[] gerpScoreArray = new float[830829764];

            int pos = -1;
            while ((temp = brg.readLine()) != null) {
                tem = KStringUtils.fastSplitdel(temp);
                pos++;
                gerpScoreArray[pos] = Float.parseFloat(tem[1]);
            }
            brg.close();

            while ((temp = brb.readLine()) != null) {
                tem = KStringUtils.fastSplitdel(temp);
                bw.write("fixedStep" + " " + "chrom=" + tem[0] + subGenome + " " + "start=" + String.valueOf(Integer.parseInt(tem[1]) + 1) + " " + "step=1");
                for (int i = Integer.parseInt(tem[1]); i < Integer.parseInt(tem[2]); i++) {
                    bw.newLine();
                    bw.write(String.valueOf(gerpScoreArray[i]));
                }
                bw.newLine();
            }

            brb.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in convert IWGSC GERP.rates to phyloP format! please try again! Come on !");
            e.printStackTrace();
        }
    }

//mafRevise___convert bed getfasta positions in msa to regular msa
//example:  
//a score=278 mismap=2.38e-06
//s traesA.chr3                100018229 74 + 750843639 atGGTCACCTTGAGGTAGGTGTCGACGGCGCGGTAGAGCGCGTCgtcggcggcgCGCGcgtgggcggGCAcagc
//s leper.chr1:3115495-3195902     80333 74 -     80407 ACGCTGACCTTGAGGTAGGTGTCGACGGCGCGGTAGAGCCCATCATCGGCGGGCCGCGCGTGGGCCGGCACGGC 
//##s leper.chr1:x1-x2  x3  x4  -   x5  .....
    //lper.chr1 realChrlength-x2+x3    x4  - readChrLength  ......   

    public void mafRevise(String mafDir, String outDir, String specie, String faiDir) {
      BufferedReader brm ;
        BufferedReader brf = YaoIOUtils.getTextReader(faiDir);
        BufferedWriter bw = YaoIOUtils.getTextWriter(outDir);
        HashMap<String, Integer> chrInfo = new HashMap();
             if (mafDir.endsWith("gz")) {
            brm = YaoIOUtils.getTextGzipReader(mafDir);
        } else {
            brm = YaoIOUtils.getTextReader(mafDir);
        }
        try {
            while ((temp = brf.readLine()) != null) {
                chrInfo.put(KStringUtils.fastSplitdel(temp)[0], Integer.valueOf(KStringUtils.fastSplitdel(temp)[1]));
            }
            brf.close();
            while ((temp = brm.readLine()) != null) {
                if (temp.startsWith("s " + specie)) {
                    String tem[] = KStringUtils.mafSplit(temp);
                    String[] info = KStringUtils.fastSplitColon(tem[1]);
                    bw.write(tem[0] + "\t");
                    bw.write(info[0] + "\t");
                    if(tem[4].equals("-")){
                    bw.write(String.valueOf(Integer.valueOf(chrInfo.get(info[0]))-Integer.valueOf(KStringUtils.fastSplitStrigula(info[1])[1]) + Integer.valueOf(tem[2])) + "\t");
                    }
                    if(tem[4].equals("+")){
                    bw.write(String.valueOf(Integer.valueOf(KStringUtils.fastSplitStrigula(info[1])[0]) + Integer.valueOf(tem[2])) + "\t");
                    }
                    bw.write(tem[3] + "\t");
                    bw.write(tem[4] + "\t");
                    Integer m = chrInfo.get(info[0]);
                    String x = String.valueOf(chrInfo.get(info[0]));
                    bw.write(String.valueOf(chrInfo.get(info[0])) + "\t");
                    bw.write(tem[6] + "\t");
                    bw.newLine();
                } else {
                    bw.write(temp);
                    bw.newLine();
                }
            }
            brm.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in revise maf! please try again! Come on !");
            e.printStackTrace();
        }
    }

    public void order(String mafDir, String outDir) {

        BufferedReader brb = YaoIOUtils.getTextReader(mafDir);
        BufferedWriter bw = YaoIOUtils.getTextWriter(outDir);
        int chr = 0;
        int chrABreakpoint[] = {471304005, 462376173, 454103970, 452555092, 453230519, 452440856, 450046986};
        int chrBBreakpoint[] = {438720154, 453218924, 448155269, 451014251, 451372872, 452077197, 453822637};
        int chrDBreakpoint[] = {452179604, 462216879, 476235359, 451004620, 451901030, 450509124, 453812268};
        int start = 0;
        int end = 0;
        String subGeno = null;
        int breakPoint[] = new int[7];
        try {
            while ((temp = brb.readLine()) != null) {
                if (temp.startsWith("chrUn")) {
                    bw.write(temp);
                    bw.newLine();
                    continue;
                }
                tem = KStringUtils.fastSplitdel(temp);
                chr = Integer.parseInt(tem[0].substring(3, 4));
                subGeno = tem[0].substring(4, 5);
                if (subGeno.equals("A")) {
                    breakPoint = chrABreakpoint;
                } else if (subGeno.equals("B")) {
                    breakPoint = chrBBreakpoint;
                } else if (subGeno.equals("D")) {
                    breakPoint = chrDBreakpoint;
                }
                if (tem[0].endsWith("part2")) {
                    start = Integer.parseInt(tem[1]) + breakPoint[chr - 1];
                    end = Integer.parseInt(tem[2]) + breakPoint[chr - 1];
                    bw.write("chr" + chr + subGeno + "\t" + start + "\t" + end);
                    for (int i = 3; i < tem.length; i++) {
                        bw.write("\t" + tem[i]);
                    }
                    bw.newLine();
                } else {
                    bw.write("chr" + chr + subGeno + "\t" + tem[1] + "\t" + tem[2]);
                    for (int i = 3; i < tem.length; i++) {
                        bw.write("\t" + tem[i]);
                    }
                    bw.newLine();
                }
            }

            brb.close();

            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in order bed! please try again! Come on !");
            e.printStackTrace();
        }
    }

    public void est(String mafDir, String outDir, String subGenome) {
        BufferedReader brb = YaoIOUtils.getTextReader(mafDir);
        BufferedWriter bw = YaoIOUtils.getTextWriter(outDir);
        String species[] = {subGenome, "secer", "hovul"};
        Byte outGroup[] = new Byte[2];
        String pos = null;
        String chr = null;
        String specie = null;
        try {
            while ((temp = brb.readLine()) != null) {
                if (temp.startsWith("a")) {

                    if (chr == null) {
                        continue;
                    }
                    bw.write(chr + "\t");
                    bw.write(pos + "\t");
                    bw.write(String.valueOf(Integer.parseInt(pos) + 1));
                    for (int i = 0; i < outGroup.length; i++) {
                        bw.write("\t");
                        if (outGroup[i] == null) {
                            bw.write("0,0,0,0");
                        } else {
                            if (outGroup[i] == 0) {
                                bw.write("1,0,0,0");
                            } else {
                                if (outGroup[i] == 1) {
                                    bw.write("0,1,0,0");
                                } else {
                                    if (outGroup[i] == 2) {
                                        bw.write("0,0,1,0");
                                    } else {
                                        if (outGroup[i] == 3) {
                                            bw.write("0,0,0,1");
                                        }
                                    }

                                }

                            }
                        }
                    }
                    bw.newLine();
                    outGroup[0] = null;
                    outGroup[1] = null;
                    //  outGroup[2] = null;
                    chr = null;
                    pos = null;
                }

                if (temp.startsWith("s")) {

                    tem = KStringUtils.mafSplit(temp);
                    if (tem[1].startsWith(subGenome)) {
                        chr = KStringUtils.fastSplitDot(tem[1])[1];
                        pos = tem[2];
                    }
                    specie = KStringUtils.fastSplitDot(tem[1])[0];
                    for (int i = 1; i < species.length; i++) {
                        if (species[i].equals(specie)) {
                            outGroup[i - 1] = KStringUtils.baseToByte(tem[6]);

                        }
                    }
                }
            }

            bw.write(chr + "\t");
            bw.write(pos + "\t");
            bw.write(String.valueOf(Integer.parseInt(pos) + 1));
            for (int i = 0; i < outGroup.length; i++) {
                bw.write("\t");
                if (outGroup[i] == null) {
                    bw.write("0,0,0,0");
                } else {
                    if (outGroup[i] == 0) {
                        bw.write("1,0,0,0");
                    } else {
                        if (outGroup[i] == 1) {
                            bw.write("0,1,0,0");
                        } else {
                            if (outGroup[i] == 2) {
                                bw.write("0,0,1,0");
                            } else {
                                if (outGroup[i] == 3) {
                                    bw.write("0,0,0,1");
                                }
                            }

                        }

                    }
                }
            }
            brb.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in order bed! please try again! Come on !");
            e.printStackTrace();
        }
    }

    public void focal(String mafDir, String outDir, String subGenome) {
        BufferedReader brb = YaoIOUtils.getTextReader(mafDir);
        BufferedWriter bw = YaoIOUtils.getTextWriter(outDir);
        //  Byte focal[] = new Byte[840000000];
        Byte focal[] = new Byte[84000];
        int pos = 0;
        int cn = 0;
        String chr = null;
        String specie = null;
        try {
            while ((temp = brb.readLine()) != null) {
                if (temp.startsWith("s")) {
                    tem = KStringUtils.mafSplit(temp);
                    if (tem[1].startsWith(subGenome)) {
                        chr = KStringUtils.fastSplitDot(tem[1])[1].substring(3, 4);
                        pos = Integer.parseInt(tem[2]);

                        for (int j = 0; j < tem[6].length(); j++) {
                            if (String.valueOf(tem[6].charAt(j)).equals("-")) {
                                continue;
                            }
                            focal[pos] = KStringUtils.baseToByte(String.valueOf(tem[6].charAt(j)));
                            pos++;
                        }
                    }
                }
            }
            for (int i = 0; i < focal.length; i++) {
                if (focal[i] == null) {
                    continue;
                }
                bw.write("chr" + chr);
                bw.write("\t");
                bw.write(String.valueOf(i));
                bw.write("\t");
                bw.write(String.valueOf(i + 1));
                bw.write("\t");

                if (focal[i] == 0) {
                    bw.write("20,0,0,0");
                } else {
                    if (focal[i] == 1) {
                        bw.write("0,20,0,0");
                    } else {
                        if (focal[i] == 2) {
                            bw.write("0,0,20,0");
                        } else {
                            if (focal[i] == 3) {
                                bw.write("0,0,0,20");
                            }
                        }
                    }

                }
                bw.newLine();
            }

            brb.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in order bed! please try again! Come on !");
            e.printStackTrace();
        }
    }

    public void dfe(String mafDir, String outDir, String subGenome) {
        BufferedReader brb = YaoIOUtils.getTextReader(mafDir);
        BufferedWriter bw = YaoIOUtils.getTextWriter(outDir);
        String species[] = {"secer", "hovul", subGenome};
        Byte outGroup[][] = new Byte[2][840000000];
        // Byte outGroup[][] = new Byte[2][84000];
        int pos = 0;
        int cn = 0;
        String chr = null;
        String specie = null;
        ArrayList<Integer> indel = new ArrayList<>();
        try {
            while ((temp = brb.readLine()) != null) {
                if (temp.startsWith("s")) {
                    tem = KStringUtils.mafSplit(temp);
                    if (tem[1].startsWith(subGenome)) {
                        indel = new ArrayList<>();
                        chr = KStringUtils.fastSplitDot(tem[1])[1].substring(3, 4);
                        pos = Integer.parseInt(tem[2]);
                        for (int j = 0; j < tem[6].length(); j++) {
                            if (String.valueOf(tem[6].charAt(j)).equals("-")) {
                                indel.add(j);
                                continue;
                            }
                        }
                    }
                    specie = KStringUtils.fastSplitDot(tem[1])[0];
                    cn = pos;
                    for (int i = 0; i < species.length - 1; i++) {
                        if (species[i].equals(specie)) {
                            for (int j = 0; j < tem[6].length(); j++) {
                                if (indel.contains(j)) {
                                    continue;
                                }
                                if (String.valueOf(tem[6].charAt(j)).equals("-")) {
                                    cn++;
                                    continue;
                                }
                                outGroup[i][cn] = KStringUtils.baseToByte(String.valueOf(tem[6].charAt(j)));
                                cn++;
                            }
                        }
                    }
                }
            }
            for (int i = 0; i < outGroup[1].length; i++) {
                if (outGroup[1][i] == null || outGroup[0][i] == null) {
                    continue;
                }
                bw.write("chr" + chr);
                bw.write("\t");
                bw.write(String.valueOf(i));
                bw.write("\t");
                bw.write(String.valueOf(i + 1));
                bw.write("\t");
                for (int j = 0; j < species.length - 1; j++) {
                    if (outGroup[j][i] == 0) {
                        bw.write("A");
                    } else {
                        if (outGroup[j][i] == 1) {
                            bw.write("C");
                        } else {
                            if (outGroup[j][i] == 2) {
                                bw.write("G");
                            } else {
                                if (outGroup[j][i] == 3) {
                                    bw.write("T");
                                }
                            }

                        }

                    }
                    if (j == 0) {
                        bw.write("\t");
                    }
                    if (j == 1) {
                        bw.newLine();
                    }
                }

            }
            brb.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in order bed! please try again! Come on !");
            e.printStackTrace();
        }
    }

    public void filterINSIGHT(String mafDir, String outDir) {
        BufferedReader brb = YaoIOUtils.getTextReader(mafDir);
        BufferedWriter bw = YaoIOUtils.getTextWriter(outDir);
        ArrayList<String> sites = new ArrayList<>();
        ArrayList<String> block = new ArrayList<>();
        try {
            temp = brb.readLine();
            bw.write(temp);
            while ((temp = brb.readLine()) != null) {
                tem = KStringUtils.fastSplitdel(temp);
                if (temp.startsWith("block")) {

                    if (block.size() == 1) {
                        bw.newLine();
                        bw.write(block.get(0));
                        for (int x = 0; x < sites.size(); x++) {
                            bw.newLine();
                            bw.write(sites.get(x));
                        }
                    }
                    if (temp.equals("block")) {
                        continue;
                    }
                    sites = new ArrayList<String>();
                    block = new ArrayList<String>();
                    if (Double.parseDouble(tem[5]) < 0.0005) {
                        block.add(temp);

                    } else {

                        brb.readLine();
                    }

                } else {
                    sites.add(temp);
                }

            }

            brb.close();

            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in order bed! please try again! Come on !");
            e.printStackTrace();
        }
    }

    public void tata(String mafDir, String outDir) {
        BufferedReader brb = YaoIOUtils.getTextReader(mafDir);
        BufferedWriter bw = YaoIOUtils.getTextWriter(outDir);
        String chr = null;
        int start = 0;
        int end = 0;
        try {
            while ((temp = brb.readLine()) != null) {
                if (temp.startsWith("Query")) {
                    tem = KStringUtils.fastSplitSpa(temp).toArray(new String[KStringUtils.fastSplitSpa(temp).size()]);
                    chr = KStringUtils.fastSplitColon(tem[1])[0].substring(1, 6);
                    start = Integer.parseInt(KStringUtils.fastSplitStrigula(KStringUtils.fastSplitColon(tem[1])[1])[0]);
                    end = Integer.parseInt(KStringUtils.fastSplitStrigula(KStringUtils.fastSplitColon(tem[1])[1])[1]);
                }
                if (temp.startsWith("TSS")) {
                    bw.write(chr + "\t");
                    bw.write(start + "\t");
                    bw.write(end + "\t");
                    bw.write(KStringUtils.mafSplit(temp)[2] + "\t");
                    bw.write(KStringUtils.mafSplit(temp)[6] + "\t");
                    bw.write(KStringUtils.mafSplit(temp)[9] + "\t");
                    bw.write(KStringUtils.mafSplit(temp)[13]);
                    bw.newLine();
                }
                if (temp.startsWith("(TATA-)")) {
                    bw.write(chr + "\t");
                    bw.write(start + "\t");
                    bw.write(end + "\t");
                    tem = KStringUtils.mafSplit(temp);
                    bw.write(KStringUtils.mafSplit(temp)[3] + "\t");
                    bw.write(KStringUtils.mafSplit(temp)[7]);
                    bw.newLine();
                }
            }
            brb.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in order bed! please try again! Come on !");
            e.printStackTrace();
        }
    }

    public void createOptions() {
        options = new Options();
        options.addOption("m", true, "maf path");
        options.addOption("f", true, "fai Path");
        options.addOption("o", true, "output Path");
        options.addOption("s", true, "specie Path");
        options.addOption("g", true, "subGenome Path");
    }

    public void retrieveParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            mafDir = line.getOptionValue("m");
            faiDir = line.getOptionValue("f");
            outDir = line.getOptionValue("o");
            specie = line.getOptionValue("s");
            subGenome = line.getOptionValue("g");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        new FormatConvertTools(args);
    }
}
