/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import com.koloboke.collect.map.hash.HashIntIntMap;
import static com.koloboke.collect.map.hash.HashIntIntMaps.newMutableMap;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import lipengKang.analysis.KStringUtils;
import utils.PStringUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author kanglipeng
 */
public class LongReadQC {

    public LongReadQC() {
        this.ctgQC();
        //this.geneRange();
    }
            String temp = null;
        String tem[] = null;

    private void geneRange() {
        String gff3Dir = "/data1/home/lipeng/database/genome/urartu/WheatTu.annotation/WheatTu.gene.gff";
        //   String gff3Dir ="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.gff3";
        //String gff3Dir="/data1/home/lipeng/database/genome/wheat/annotation/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gff3";
        BufferedReader br;
        if (gff3Dir.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(gff3Dir);
        } else {
            br = YaoIOUtils.getTextReader(gff3Dir);
        }

        int geneStartPos = 0;
        int geneEndPos = 0;
        int geneIndex = 0;
        //hashmap: chrNum------>sort.geneRanges
        HashIntIntMap[] geneRangesMap = new HashIntIntMap[7];
        for (int i = 0; i <= 6; i++) {
            geneRangesMap[i] = newMutableMap();
        }

        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("TuUngrouped")) {
                    continue;
                }

                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
//                    StringBuilder chr = new StringBuilder();`
//                    chr.append(String.valueOf(i)).append("A");
//                    int chr = Integer.valueOf(String.valueOf(tem[0].charAt(3)));

                if (tem[2].equals("gene")) {
                    geneIndex++;
                    int chr = Integer.valueOf(String.valueOf(tem[0].charAt(2)));
                    geneStartPos = Integer.parseInt(tem[3]);
                    geneEndPos = Integer.parseInt(tem[4]);
                    for (int j = geneStartPos; j <= geneEndPos; j++) {
                        geneRangesMap[chr - 1].put(j, geneIndex);
                    }
                }
            }

            br.close();

            System.out.println("indexing genes in urartu genome. ok!");
        } catch (Exception e) {
            System.out.println("Error in indexing genes in urartu genome");
            e.printStackTrace();
        }

        String pafFile = "/data1/home/lipeng/result/PacBioAnalysis/urartuQC1/readschain.kbmap";
        //    String pafFile ="/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.paf";
        BufferedReader bs;
        if (pafFile.endsWith("gz")) {
            bs = YaoIOUtils.getTextGzipReader(pafFile);
        } else {
            bs = YaoIOUtils.getTextReader(pafFile);
        }
        int mappedBlockStart = 0;
        int mappedBlockEnd = 0;
        long contigsSumLength = 0;
        int ctgIndex = 0;
        HashIntIntMap[] ctgRangesMap = new HashIntIntMap[7];
        for (int i = 0; i <= 6; i++) {
            ctgRangesMap[i] = newMutableMap();
        }
        try {
            while ((temp = bs.readLine()) != null) {

                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (tem[5].startsWith("TuU")) {
                    continue;
                }
                ctgIndex++;
                mappedBlockStart = Integer.parseInt(tem[8]);
                mappedBlockEnd = Integer.parseInt(tem[9]);
                contigsSumLength = contigsSumLength + mappedBlockEnd - mappedBlockStart;
                int chr = Integer.valueOf(String.valueOf(tem[5].charAt(2)));
                for (int j = mappedBlockStart; j <= mappedBlockEnd; j++) {
                    ctgRangesMap[chr - 1].put(j, ctgIndex);
// 
                }
            }
            System.out.println("indexing wtdbg2-contigs. ok!");
            bs.close();
            //create mapped contigs pos lists for each chr 
//            for (int i = 0; i <= 6; i++) {
//               for(int j: ctgRangesMap[i].keySet()){
//               geneRangesMap[i].keySet().
//               
//               }
//      
//          System.out.println(" wtdbg2-contigs index map constructing. ok!");

        } catch (Exception e) {
            System.out.println("Error in indexing genes in urartu genome");
            e.printStackTrace();
        }
        //  return genePosList;
    }
    
    
    //------------------------wtdbg2+wtpoa-->ctg.fa--------------------------------------------
    //ctg.fa convert by grep '>ctg' --->ctg.list 
    // ctg.list format example----->  >ctg1 len=268589 /n >ctg2 len=267931


    public void ctgQC(){
     String ctgListDir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/homo.ctg.list";
    //String ctgListDir = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/gerp1.txt";
        BufferedReader br;
        if (ctgListDir.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(ctgListDir);
        } else {
            br = YaoIOUtils.getTextReader(ctgListDir);
        }
        ArrayList <Integer> ctgList=new ArrayList();
         try {
            while ((temp = br.readLine()) != null) {
                  List<String> tList = KStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                ctgList.add(Integer.valueOf(tem[1].substring(4)));
                
                
        
         }
            br.close();
         Collections.sort(ctgList);
         Long sumCtg=0L;
         for(int i:ctgList){
             sumCtg=sumCtg+i;
         }
         Long sum1 =0L;
         Long sum2 =0L;
         for(int j=0;j<ctgList.size()-1;j++){
             sum1=sum1+ctgList.get(j);
             sum2=sum1+ctgList.get(j+1);
        if (sum1<sumCtg/2 &&sum2>=sumCtg/2){
       System.out.println("N50 "+ ctgList.get(j+1));
        }
        if (sum1==sumCtg/2 ){
       System.out.println("N50 "+ ctgList.get(j));
        }
        
         }
         System.out.println("MIN "+ctgList.get(0));
          System.out.println("MAX "+ctgList.get(ctgList.size()-1));
           System.out.println("AVG "+sumCtg/ctgList.size());
           
           System.out.println("output "+ctgList.size()+" ctgs");
           System.out.println("TOT "+sumCtg);
           
         
         }catch (Exception e) {
            System.out.println("Error in indexing genes in urartu genome");
            e.printStackTrace();
        }
    
    }
    
    

    public static void main(String[] args) {

        new LongReadQC();
    }

}
