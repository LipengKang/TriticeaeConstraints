/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.List;
import utils.PStringUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author kanglipeng
 */
public class LongReadQC {
    public LongReadQC(){}
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
        String temp=null;
        int geneStartPos = 0;
        int geneEndPos = 0;
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
            }
            
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
    
    
    
    
}
