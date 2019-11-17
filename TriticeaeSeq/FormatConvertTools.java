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
import lipengKang.analysis.KStringUtils;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
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
        // this.debugGERP();
        //   this.mafToFa();
       // this.mafTophy();
    this.gerpReformat();
        //this.mafRevise(mafDir,outDir,specie,faiDir);
    }
//---------------------GTF to siftGTF--------------------------------------------------
    //convert gtf from IWGSC to sift4g input gtf
    //sift-gtf requires 'gene_biotype "protein_coding"' in ninth column of normal gtf

    public void debugGERP() {

        String GTF = "/data1/home/lipeng/result/GERP/axt/A/gerp/wheatA.maf";
        BufferedReader br;

        //info point to ninth column of GTF
        if (GTF.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(GTF);

        } else {
            br = YaoIOUtils.getTextReader(GTF);

        }

        try {

            int gerpPos = 0;
            //extract info of each transcript   
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("s wheatA")) {
                    //  if (temp.startsWith("s Tu1")) {
                    List<String> fList = KStringUtils.fastSplit(temp);
                    List<String> fListNew = new ArrayList<>();
                    for (int i = 0; i < fList.size(); i++) {
                        if (fList.get(i) != null && !fList.get(i).equals("")) {
                            fListNew.add(fList.get(i));
                        }
                    }
                    tem = fListNew.toArray(new String[fListNew.size()]);

                    gerpPos = gerpPos + tem[6].replaceAll("N", "").replaceAll("n", "").replaceAll("-", "").length();
                    if (tem[2].equals("508641235")) {
                        break;
                    }
                }

            }

            br.close();
            System.out.println(gerpPos);

        } catch (Exception e) {
            System.out.println("Error in convert IWGSC GTF to sift4g GTF! please try again! Come on !");
            e.printStackTrace();
        }
        System.out.println("Congratulations! IWGSC GTF converts to sift4g GTF sucessfully!");
    }

    //mafToFa v1.1 : convert maf format to fa format. pay attention! this convertor only used to bring up fa format for GERP++
    public void mafToFa() {
        String maf = "/data2/lipeng/msa/multiz/a.4d.maf";
        String out = "/data2/lipeng/msa/multiz/a.4d.fa";
        BufferedReader br = YaoIOUtils.getTextReader(maf);
        BufferedWriter bw = YaoIOUtils.getTextWriter(out);

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

    //maf2phy
    public void mafTophy() {
        String maf = "/Users/kanglipeng/Desktop/test.maf";
        String out = "/Users/kanglipeng/Desktop/test.phy";
        BufferedReader br = YaoIOUtils.getTextReader(maf);
        BufferedWriter bw = YaoIOUtils.getTextWriter(out);
        int specieNum = 0;
        String fourDbase = null;
        String[] species = {"traes", "aetau", "trura", "hospo", "hovul", "brdis", "phedu", "leper", "orbra", "orpun", "ormer", "orlon", "orglu", "orbar", "orgla", "orniv", "orind", "orruf", "orjap", "ertef", "ercur", "ortho", "sobic", "zemay", "seita", "eccru", "pahal", "pavir", "ancom"};
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
            System.out.println("Error in convert IWGSC GTF to sift4g GTF! please try again! Come on !");
            e.printStackTrace();
        }
        System.out.println("Congratulations! IWGSC GTF converts to sift4g GTF sucessfully!");
    }

    //gerpNorm: convert GERP.rates to phylop standard output format
    //example:
// GERP.rates                       phylop format(1-based)
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
            float[] gerpScoreArray=new float[830829764];

            int pos = -1;
            while ((temp = brg.readLine()) != null) {
                tem = KStringUtils.fastSplitdel(temp);
                pos++;
                gerpScoreArray[pos]= Float.parseFloat(tem[1]);
            }
              brg.close();
      
            while ((temp = brb.readLine()) != null) {
                tem = KStringUtils.fastSplitdel(temp);
                bw.write("fixedStep"+" "+"chrom="+tem[0]+subGenome+" "+"start="+String.valueOf(Integer.parseInt(tem[1])+1)+" "+"step=1");
                for(int i=Integer.parseInt(tem[1]);i<Integer.parseInt(tem[2]);i++){
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

    public void mafRevise(String mafDir, String outDir, String specie, String faiDir) {

        BufferedReader brm = YaoIOUtils.getTextReader(mafDir);
        BufferedReader brf = YaoIOUtils.getTextReader(faiDir);
        BufferedWriter bw = YaoIOUtils.getTextWriter(outDir);
        HashMap<String, Integer> chrInfo = new HashMap();

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
                    bw.write(String.valueOf(Integer.valueOf(KStringUtils.fastSplitStrigula(info[1])[0]) + Integer.valueOf(tem[2])) + "\t");
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
