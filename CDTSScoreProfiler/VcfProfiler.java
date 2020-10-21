/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CDTSScoreProfiler;

import com.koloboke.collect.map.hash.HashIntFloatMap;
import com.koloboke.collect.map.hash.HashIntFloatMaps;
import format.dna.FastaByte;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import TriticeaeSeq.KStringUtils;

/**
 *
 * @author lipeng
 * extract position and AF of all training datasets(VCFs)
 */
public class VcfProfiler {
     HashIntFloatMap[] AFScoreMaps = null;
//    public VcfProfiler(String vcfDirs,FastaByte f){
//    this.calculateAF(vcfDirs,f);
//    }

    public static HashIntFloatMap[] calculateAF(String vcfDirS,FastaByte f) {
        System.out.println("vcf profiling...");
        String temp = null;
        String[] tem = null;
        String[] infoTem = null;
      //  int DP = 0;
        String ADString[] = null;
        float AFScore = 0f;
        int chr = 0;
       
        //maps defined by chrs
        HashIntFloatMap[] AFScoreMaps = new HashIntFloatMap[f.getSeqNumber()];
        
        for (int i = 0; i < f.getSeqNumber(); i++) {
            AFScoreMaps[i] = HashIntFloatMaps.newMutableMap(f.getSeqLength(i));
        }
       
        File file = new File(vcfDirS);

File[] files = file.listFiles(); 
 List<String>list=new ArrayList ();
 System.out.println("list vcf files...");
   for (int i=0;i<files.length;i++) {                    
System.out.println("vcf file "+String.valueOf(i+1)+" : " + files[i].getAbsolutePath());
list.add(files[i].getAbsolutePath());
}
    String vcfs[]=new String[list.size()];
    list.toArray(vcfs);
    System.out.println("calculate AF from vcf files......");

        for(String vcf:vcfs){
        System.out.println("dealing with : " + vcf+" ......")
;        BufferedReader br;
        if (!vcf.endsWith(".gz")) {
            br = utils.IOUtils.getTextReader(vcf);
        } else {
            br = utils.IOUtils.getTextGzipReader(vcf);
        }

        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    continue;
                }
                //split for ADs and DP
//??????                List<String> fList = KStringUtils.fastSplitdel(temp);
  //?? ??            tem = fList.toArray(new String[fList.size()]);
                List<String> sList = KStringUtils.fastSplitSemicolon(tem[7]);
                infoTem = sList.toArray(new String[sList.size()]);
                //AFScore calculate
             //   DP = Integer.valueOf(infoTem[0].substring(3));
             //    ADString = KStringUtils.fastSplitComma(infoTem[1].substring(3));
               
//                int AD[] = new int[ADString.length];
//                for (int i = 0; i < ADString.length; i++) {
//                    AD[i] = Integer.valueOf(ADString[i]).intValue();

//extract MAF
AFScore = Float.parseFloat(infoTem[6].substring(4));

            //    }
              //  Arrays.sort(AD);
           //      DecimalFormat df=new DecimalFormat("0.000");
          //      AFScore=Float.parseFloat(df.format((float)AD[AD.length - 1] / DP));
                //chrmosome name convert like 1A->1,2A->2,1B->8,1D->15
//                if (tem[0].endsWith("A")) {
//                    chr = -1+Integer.valueOf(tem[0].substring(0, 1));
//                }
//                if (tem[0].endsWith("B")) {
//                    chr = 6 + Integer.valueOf(tem[0].substring(0, 1));
//                }
//                if (tem[0].endsWith("D")) {
//                    chr = 13 + Integer.valueOf(tem[0].substring(0, 1));
//                }
// chromosme name convert from lulab wheat fasta format.
//chr1 -->AFScoreMaps[0]. 
chr=Integer.valueOf(tem[0])-1;
if(chr<f.getSeqNumber()){
                //vcf is 1based, but AFScoreMas is 0-based
               AFScoreMaps[chr].put(Integer.valueOf(tem[1]).intValue()-1, AFScore);}
            }
            br.close();
        } catch (Exception e) {
            System.out.println("Errors in vcf profiling! please try again! Come on !");
                e.printStackTrace();
        }}
        System.out.println("AFScores calculation finished");
        return AFScoreMaps;
    }
}
