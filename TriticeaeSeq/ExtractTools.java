/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import format.dna.BaseEncoder;
import format.dna.FastaByte;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import lipengKang.analysis.KStringUtils;
import static lipengKang.analysis.RandomArray.randomArray;
import utils.IOUtils;
import utils.PStringUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

public class ExtractTools {

    public ExtractTools() {
          this.genomeToReads();
//        this.mafToBed();
    }

    //----------------cut genome to n x accurate 100% reads----------------------
    //ensure consistency of triticeae alignment proccess
    public void genomeToReads() {
//       hovul.chr1	558535432
//hovul.chr2	768075024
//hovul.chr3	699711114
//hovul.chr4	647060158
//hovul.chr5	670030160
//hovul.chr6	583380513
//hovul.chr7	657224000
//hovul.chrUn	249774706
int[]chrIndex={558535432,768075024,699711114,647060158,670030160,583380513,657224000,249774706};
        try {
          BufferedWriter bw = YaoIOUtils.getTextWriter("/Users/kanglipeng/Desktop/4.bed");

            for(int i=0;i<=7;i++){
               int chr=i+1;
               bw.write("hovul.chr"+chr+"\t"+0+"\t"+14999);
            bw.newLine();
            int start =15000;int end=34999;
            bw.write("hovul.chr"+chr+"\t"+start+"\t"+end);
            bw.newLine();
            while(start<chrIndex[i]-20000&&end<=chrIndex[i]-20000){
            start=end+1;end=start+19999;
             bw.write("hovul.chr"+chr+"\t"+start+"\t"+end);
                  bw.newLine();}
            start=end+1;end=start+19999;
            if(start>chrIndex[i])continue;
            if(start<=chrIndex[i]&&end>chrIndex[i]){ bw.write("hovul.chr"+chr+"\t"+start+"\t"+chrIndex[i]);
            if(i!=7){
                  bw.newLine();}}
            }
bw.flush();
bw.close();
        } catch (Exception e) {
            System.out.println("Errors in vcf profiling! please try again! Come on !");
            e.printStackTrace();
        }

    }

    public void mafToBed() {
         String fa = "/data1/home/lipeng/database/assembly/barley5x/psvil.unassembled.fasta";
          String outfileS = "/data1/home/lipeng/database/assembly/barley5x/psvil.flt.unassembled.fasta";
//        String fa = "/Users/kanglipeng/Desktop/test.fa";
//        String outfileS = "/Users/kanglipeng/Desktop/test.txt";
        try {
            String temp = null;
            String tem[] = null;
            BufferedReader br = YaoIOUtils.getTextReader(fa);
            BufferedWriter bw = YaoIOUtils.getTextWriter(outfileS);
            temp = br.readLine();
            while (temp != null && temp.startsWith(">")) {
                List<String> tList = KStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (Integer.valueOf(tem[2].substring(6)) >= 2) {
                    bw.write(temp);
                    bw.newLine();
                    while ((temp = br.readLine()) != null && !temp.startsWith(">")) {
                        bw.write(temp);
                        bw.newLine();
                    }
                } else {
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith(">")) {
                            break;
                        }
                    }
                }
            }

            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            System.out.println("Errors in vcf profiling! please try again! Come on !");
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        new ExtractTools();
    }

}
