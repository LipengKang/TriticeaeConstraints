/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashIntByteMap;
import com.koloboke.collect.map.hash.HashIntByteMaps;
import format.dna.BaseEncoder;
import format.dna.FastaByte;
import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import lipengKang.analysis.KStringUtils;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import utils.IOUtils;

/**
 *
 * @author kanglipeng
 */
public class INSIGHT {

    Options options = new Options();
    String temp = null;
    String tem[] = null;
    String lambdaFile = null;
    String probFile = null;
    String neutralFile = null;
    String outFile = null;
    String faFile = null;
    String vcfFile = null;

    public INSIGHT(String[] args) {
        this.createOptions();
        this.retrieveParameters(args);
       // this.formatINS();
        this.filterIns();
    }

//requirements:1.lambda two column file "blockID    lambda" 2.
    // public void formatINS(String lambdaDir, String prequelFile, String snpFile, String neutralFile, String fa,String outDir) {
    public void formatINS() {

        BufferedReader brl = IOUtils.getTextReader(lambdaFile);
        BufferedReader brp = IOUtils.getTextReader(probFile);
        // BufferedReader brs = IOUtils.getTextGzipReader(vcfFile);
        BufferedReader brs = IOUtils.getTextReader(vcfFile);
        BufferedReader brn = IOUtils.getTextReader(neutralFile);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        FastaByte f = new FastaByte(faFile);
        HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIByteMap();
        //only first seq in fa  ACGT-->(0123)
        byte[] bArray = f.getSeq(0).getBytes();
        for (int j = 0; j < bArray.length; j++) {
            bArray[j] = ascIIByteMap.get(bArray[j]);
        }
        try {
            //process lambda
            //format:block_NUM  lambda
            String lambda = null;
            HashMap<Integer, String> lambdaMap = new HashMap<>();
            while ((temp = brl.readLine()) != null) {
                tem = KStringUtils.fastSplitdel(temp);
                if(Double.parseDouble(tem[1])>0.00148569993756236178d){continue;}
                lambda = String.format("%10.9e", Double.parseDouble(tem[1]));
                lambdaMap.put(Integer.parseInt(tem[0]), lambda);
            }
            brl.close();
            //initiate neurtal pos---> ref       M:alt=0 P:alt=1234(ACGT)
            //   HashIntByteMap neutralSitesInfo = HashIntByteMaps.newMutableMap();
            int[] neutralIndex = new int[850000000];
            for (int i = 0; i < neutralIndex.length; i++) {
                neutralIndex[i] = 4;
            }
            while ((temp = brn.readLine()) != null) {
                tem = KStringUtils.fastSplitdel(temp);
                for (int i = Integer.parseInt(tem[1]); i < Integer.parseInt(tem[2]); i++) {
                    //       neutralSitesInfo.put(i, bArray[i]);
                    neutralIndex[i] = bArray[i];
                }
            }
            brn.close();
            System.out.println("lambda map established");
            //add  snps info ---> L:0 H:1   majMap-->major allele
            HashMap<Integer, String> snpInfo = new HashMap<>();
            HashIntByteMap majMap = HashIntByteMaps.newMutableMap();
            HashIntByteMap minMap = HashIntByteMaps.newMutableMap();
            float MAF = 0f;
            while ((temp = brs.readLine()) != null) {
                if (temp.startsWith("#")) {
                    continue;
                }
                tem = KStringUtils.fastSplitdel(temp);
                // AAAF_ABD>0.5 shif alt to maj
                if (tem[6] == "NaN") {
                    continue;
                }
                MAF = Float.parseFloat(tem[6]);
                //  if (neutralSitesInfo.keySet().contains(Integer.parseInt(tem[1]) - 1)) {
                if (neutralIndex[Integer.parseInt(tem[1])] != 4) {
                    if (MAF > 0.5f) {
                        MAF = 1 - MAF;
                        if (MAF <= 0.15f) {
                            snpInfo.put(Integer.parseInt(tem[1]), "L");
                            majMap.put(Integer.parseInt(tem[1]), KStringUtils.baseToByte(tem[4]));
                            minMap.put(Integer.parseInt(tem[1]), KStringUtils.baseToByte(tem[3]));
                        } else {
                            snpInfo.put(Integer.parseInt(tem[1]), "H");
                            majMap.put(Integer.parseInt(tem[1]), KStringUtils.baseToByte(tem[4]));
                            minMap.put(Integer.parseInt(tem[1]), KStringUtils.baseToByte(tem[3]));
                        }
                    } else {
                        if (MAF <= 0.15f) {
                            snpInfo.put(Integer.parseInt(tem[1]), "L");
                            majMap.put(Integer.parseInt(tem[1]), KStringUtils.baseToByte(tem[3]));
                            minMap.put(Integer.parseInt(tem[1]), KStringUtils.baseToByte(tem[4]));
                        } else {
                            snpInfo.put(Integer.parseInt(tem[1]), "H");
                            majMap.put(Integer.parseInt(tem[1]), KStringUtils.baseToByte(tem[3]));
                            minMap.put(Integer.parseInt(tem[1]), KStringUtils.baseToByte(tem[4]));
                        }
                    }
                }
            }
            brs.close();
            List<Integer> snpList = new ArrayList<>(snpInfo.keySet());
            Collections.sort(snpList);
            List<Integer> blockList = new ArrayList<>(lambdaMap.keySet());
            Collections.sort(blockList);
            System.out.println("snp majProb & minProb map established");
            // prequel 

            int pos = 0;
            HashMap<Integer, String> probMajMap = new HashMap<>();
            HashMap<Integer, String> probMinMap = new HashMap<>();
            brp.readLine();
            String majProb = null;
            String minProb = null;
            while ((temp = brp.readLine()) != null) {
                if (temp.startsWith("-")) {
                    //  neutralSitesInfo.remove(pos);
                    neutralIndex[pos] = 4;
                    pos++;
                    continue;
                }
                tem = KStringUtils.fastSplitdel(temp);
                if (tem.length != 4) {
                    neutralIndex[pos] = 4;
                    continue;
                }
                // if (neutralSitesInfo.keySet().contains(pos)) {
                if (neutralIndex[pos] != 4) {
                    // if (snpInfo.keySet().contains(pos)) {
                    if (Collections.binarySearch(snpList, pos) >= 0) {
                        //P
                        majProb = String.format("%6.5e", Double.parseDouble(tem[majMap.get(pos)]));
                        minProb = String.format("%6.5e", Double.parseDouble(tem[minMap.get(pos)]));
                        probMajMap.put(pos, majProb);
                        probMinMap.put(pos, minProb);
                    } else {
                        //M
                        majProb = String.format("%6.5e", Double.parseDouble(tem[bArray[pos]]));
                        probMajMap.put(pos, majProb);
                    }
                }
                pos++;
            }
            brs.close();
            System.out.println("ancestralProb map established");
            //theta calculation
            double monoAlleles = 0.0;
            double an = 0.0;
            double polyAlleles = 0.0;
            String theta = null;
            HashMap<Integer, String> thetaMap = new HashMap<>();
            //419 hexaploid wheat population
            for (double i = 1; i <= 418; i++) {
                an = an + 1 / i;
            }
            for (int i : blockList) {
                monoAlleles = 0.0;
                polyAlleles = 0.0;
                for (int j = i * 5000; j <= (i + 1) * 5000; j++) {
                    // if (neutralSitesInfo.keySet().contains(j)) {
                    if (neutralIndex[j] != 4) {
                        //   if (snpInfo.keySet().contains(j)) {
                        if (Collections.binarySearch(snpList, j) >= 0) {
                            polyAlleles++;
                        }
                        monoAlleles++;
                    }
                }
                if (polyAlleles == 0.0) {
                    thetaMap.put(i, "NaN");
                } else if (monoAlleles == 0.0) {
                    thetaMap.put(i, "NaN");
                } else {
                    theta = String.format("%10.9e", polyAlleles / monoAlleles / an);
                    thetaMap.put(i, theta);
                }

            }
            System.out.println("theta map established");
            System.out.println("start writing, please wait!");
            bw.write("samples" + "\t" + "838");
            bw.newLine();

            for (int i : blockList) {
                int sitesNum = 0;
                //filter 5kb blocks with netural sites <100bp (filter 15kb blocks with <100bp for estimating lambda)
                for (int j = i * 5000; j < i * 5000 + 5000; j++) {
                    if (neutralIndex[j] != 4) {
                        sitesNum++;
                    }
                }
                if (sitesNum < 200) {
                    continue;
                }
                //remove non ploly block
                if (thetaMap.get(i) == "NaN") {
                    continue;
                }
                bw.write("block" + "\t" + "chr" + f.getName(0) + ":" + i * 5000 + "-" + (i + 1) * 5000 + "\t" + "theta" + "\t" + thetaMap.get(i) + "\t" + "lambda" + "\t" + lambdaMap.get(i));
                bw.newLine();
                for (int j = i * 5000; j < i * 5000 + 5000; j++) {
                    if (neutralIndex[j] != 4) {
                        // if (snpInfo.keySet().contains(j)) {
                        if (Collections.binarySearch(snpList, j) >= 0) {
                            bw.write("site" + "\t" + "chr" + f.getName(0) + ":" + j + "\t" + snpInfo.get(j) + "\t" + probMajMap.get(j) + "\t" + probMinMap.get(j));
                            bw.newLine();
                        } else {
                            bw.write("site" + "\t" + "chr" + f.getName(0) + ":" + j + "\t" + "M" + "\t" + probMajMap.get(j));
                            bw.newLine();
                        }
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
        public void filterIns() {

        BufferedReader br = IOUtils.getTextReader(faFile);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
int blockStatus=0;
Double scale=Double.parseDouble(probFile);
//filter rule parameter =lambda/10-lambda*10
Double lambda=Double.parseDouble(lambdaFile);
        try {
            while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("sample")) {
                    bw.write(temp);bw.newLine();}
                if (temp.startsWith("block")) {
                    tem = KStringUtils.fastSplitdel(temp);
                      blockStatus=0;
                   if(Double.parseDouble(tem[5])>lambda/scale && Double.parseDouble(tem[5])<lambda*scale){
                   bw.write(temp);
                       bw.newLine();
                   blockStatus=1;
                   } 
                } 
                if (temp.startsWith("site")&&blockStatus==1) {
                    bw.write(temp);
                    bw.newLine();
                  } 
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Congratulations!");
    }

    public void createOptions() {
        options = new Options();
        options.addOption("p", true, "prequel");
        options.addOption("l", true, "lambda");
        options.addOption("o", true, "output Path");
        options.addOption("b", true, "neutral bed");
        options.addOption("f", true, "fa");
        options.addOption("v", true, "vcf");
    }

    public void retrieveParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            probFile = line.getOptionValue("p");
            lambdaFile = line.getOptionValue("l");
            faFile = line.getOptionValue("f");
            outFile = line.getOptionValue("o");
            neutralFile = line.getOptionValue("b");
            vcfFile = line.getOptionValue("v");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        new INSIGHT(args);
    }
}
