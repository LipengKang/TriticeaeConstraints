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
String Genomes[]={"A","B","D"};
//String TE[]={"DHH","DTA","DTC","DTH","DTM","DTT","DXX","RIX","RLC","RLG","RLX"};
String TE[]={"DHH","DTA","DTC","DTH","DTM","DTT"};
 int chrA[]={471304005,462376173,454103970,452555092,453230519,452440856,450046986};
        int chrB[]={438720154,453218924,448155269,451014251,451372872,452077197,453822637};
        int chrD[]={452179604,462216879,476235359,451004620,451901030,450509124,453812268};
        for (int i =0;i<=6;i++){
  
//      System.out.print("awk '{if($1==\""+"chr"+(i*6+1)+"\")print \""+(i+1)+"A"+"\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5;else print $0}' |");      
//   System.out.print("awk '{if($1==\""+"chr"+(i*6+2)+"\")print \""+(i+1)+"A\"\"\\t\"$2+"+chrA[i]+"\"\\t\"$3+"+chrA[i]+"\"\\t\"$4\"\\t\"$5;else print $0}'|");
//    System.out.print("awk '{if($1==\""+"chr"+(i*6+3)+"\")print \""+(i+1)+"B"+"\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5; else print $0}'|");      
//   System.out.print("awk '{if($1==\""+"chr"+(i*6+4)+"\")print \""+(i+1)+"B\"\"\\t\"$2+"+chrB[i]+"\"\\t\"$3+"+chrB[i]+"\"\\t\"$4\"\\t\"$5;else print $0}'|");
//    System.out.print("awk '{if($1==\""+"chr"+(i*6+5)+"\")print \""+(i+1)+"D"+"\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5;else print $0}'|");      
//   System.out.print("awk '{if($1==\""+"chr"+(i*6+6)+"\")print \""+(i+1)+"D\"\"\\t\"$2+"+chrD[i]+"\"\\t\"$3+"+chrD[i]+"\"\\t\"$4\"\\t\"$5;else print $0}'|");
//
System.out.println("grep 'chr"+(i+1)+"A' cds|awk '{if($2<"+chrA[i]+")print $0}' >chr"+(i*6+1)+".dele.bed");
System.out.println("grep 'chr"+(i+1)+"A' cds|awk '{if($2>="+chrA[i]+")print $1\"\\t\"$2-"+chrA[i]+"\"\\t\"$3-"+chrA[i]+"\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8\"\\t\"$9\"\\t\"$10\"\\t\"$11\"\\t\"$12\"\\t\"$13}' >chr"+(i*6+2)+".dele.bed");
System.out.println("grep 'chr"+(i+1)+"B' cds|awk '{if($2<"+chrB[i]+")print $0}' >chr"+(i*6+3)+".dele.bed");
System.out.println("grep 'chr"+(i+1)+"B' cds|awk '{if($2>="+chrB[i]+")print $1\"\\t\"$2-"+chrB[i]+"\"\\t\"$3-"+chrB[i]+"\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8\"\\t\"$9\"\\t\"$10\"\\t\"$11\"\\t\"$12\"\\t\"$13}' >chr"+(i*6+4)+".dele.bed");
System.out.println("grep 'chr"+(i+1)+"D' cds|awk '{if($2<"+chrD[i]+")print $0}' >chr"+(i*6+5)+".dele.bed");
System.out.println("grep 'chr"+(i+1)+"D' cds|awk '{if($2>="+chrD[i]+")print $1\"\\t\"$2-"+chrD[i]+"\"\\t\"$3-"+chrD[i]+"\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8\"\\t\"$9\"\\t\"$10\"\\t\"$11\"\\t\"$12\"\\t\"$13}' >chr"+(i*6+6)+".dele.bed");





//System.out.println("cat "+(i*5+1)+".sh "+(i*5+2)+".sh "+(i*5+3)+".sh "+(i*5+4)+".sh "+(i*5+5)+".sh >"+i+"."+"dfe"+".shs");
//for (String i :Genomes){ 
//    for (String j:TE){
   //System.out.println("grep '="+j+"' iwgsc_refseqv1.0_TransposableElements_2017Mar13.bed |awk '{if(substr($1,5,5)==\""+i+"\") print $0}' >"+i+"."+j+".bed");     
        
   // System.out.println("sed 's/A.DTH/"+i+"."+j+"/g' ffd.sh| sed 's/A.gene/"+i+".gene/g' >"+i+"."+j+".sh");
    //System.out.println("sh "+i+"."+j+".sh");
    // System.out.println("cat head >>"+i+"."+j+".dist");
    //System.out.println("paste x sort.txt >"+i+"."+j+".upstream");
        //System.out.println("awk '{print \""+i+"\"\"\\t\"\""+j+"\"\"\\t\"$2\"\\t\"$1}' "+i+"."+j+".dis >"+i+"."+j+".downstream");
//        System.out.println("bedtools flank -g wheat.genome -i "+i+"."+j+".bed -b 500 >"+i+"."+j+"_500bp.bed");
//        System.out.println("bedtools intersect -a "+i+".nc.dele -b "+i+"."+j+"_500bp.bed >"+i+"."+j+"_500bp.ncdele");
//        System.out.println("bedtools flank -g wheat.genome -i "+i+"."+j+".bed -b 1000 >"+i+"."+j+"_1000bp.bed");
//        System.out.println("bedtools intersect -a "+i+".nc.dele -b "+i+"."+j+"_1000bp.bed >"+i+"."+j+"_1000bp.ncdele");
//        System.out.println("bedtools flank -g wheat.genome -i "+i+"."+j+".bed -b 3000 >"+i+"."+j+"_3000bp.bed");
//        System.out.println("bedtools intersect -a "+i+".nc.dele -b "+i+"."+j+"_3000bp.bed >"+i+"."+j+"_3000bp.ncdele");
       //System.out.println("sortBed -i "+i+"."+j+".bed |"+"bedtools merge -i - |awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' >> class2.length");

    }
//    System.out.println("wc -l *_500bp.ncdele >500bp.ncdele");
//    System.out.println("wc -l *_1000bp.ncdele >1000bp.ncdele");
//     System.out.println("wc -l *_3000bp.ncdele >3000bp.ncdele");
//        String species="psjun";
//    System.out.println("nohup last-train -P50 --revsym --matsym --gapsym -E0.05 -C2 /data1/home/lipeng/database/lib/lastLib/NEAR/traesDb/traesDb ../../ccs/"+species+"Asse/"+species+".all.fa >../mat/traes."+species+".mat &");    
//       for (String i :Genomes){ 
//System.out.println("nohup parallel-fasta -j 10 --compress \"lastal -P6 -C2 -p /data2/lipeng/msa/mat/traes."+species+".mat /data1/home/lipeng/database/lib/lastLib/NEAR/traes"+i+"Db/traes"+i+"Db | last-split -m1 -fMAF\" < ../../ccs/"+species+"Asse/"+species+".all.fa > /data2/lipeng/msa/"+species+"/traes"+i+"."+species+".1split.maf 2>"+i+".out &");
      
    
//String annotation[]={"5UTR","intron","3UTR","core_proximal","upstream_proximal","downstream","distal","nonsyn","syn"};
// for (String i :Genomes){ 
//      for (String j :annotation){ 
//System.out.println("awk '{if($6!=$11)print \""+j+"\"\"\\t\"(1-$8);else print \""+j+"\"\"\\t\"$8}' "+i+"."+j+".daf >>"+i+".dele.daf");
//System.out.println("bedtools intersect -b "+i+"."+j+".dele -a "+i+".anc >"+i+"."+j+".daf");

      
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
