/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import lipengKang.analysis.KStringUtils;
import utils.PStringUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author kanglipeng
 */
//this package contain tools for format convert
public class FormatConvertTools {
         String temp=null;
        String[] tem=null;
        

    public FormatConvertTools() {
        this.siftGTF();
    }
//---------------------GTF to siftGTF--------------------------------------------------
    //convert gtf from IWGSC to sift4g input gtf
    //sift-gtf requires 'gene_biotype "protein_coding"' in ninth column of normal gtf
    public void siftGTF() {
        String siftGTF = "/data1/home/lipeng/sift.gtf";
        String GTF = "/data1/home/lipeng/database/genome/wheat/annotation/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gtf";
        BufferedReader br;
        BufferedWriter bw;
        //info point to ninth column of GTF
    
        if (GTF.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(GTF);
            bw = YaoIOUtils.getTextGzipWriter(siftGTF);
        } else {
            br = YaoIOUtils.getTextReader(GTF);
            bw = YaoIOUtils.getTextWriter(siftGTF);
        }
        String strand=null;
     ArrayList<Integer> cdsPosList=new ArrayList<>();
        String trsctID=null;
        Set<String>trsctSet=new HashSet<>();
         ArrayList<String>trsctList=new ArrayList<>();
        ArrayList<String>tmp=new ArrayList<>();
         ArrayList<String>trsctInfoList=new ArrayList<>();
       try {   
         String[]tems=null;
        //extract info of each transcript   
            while ((temp = br.readLine()) != null) {
                    
                    
                    List<String> tList = PStringUtils.fastSplit(temp);
                   tem = tList.toArray(new String[tList.size()]);
                   trsctID=KStringUtils.fastSplitSemicolon(tem[8]).get(0);
                   trsctSet.add(trsctID);
                   trsctList.add(trsctID);
                   tmp.add(temp);
                   
                   //transcript number>1, start parse first trancript
                   if(trsctSet.size()>1){
                   tmp.remove(tmp.size()-1);
                   trsctInfoList=tmp;
                       StringBuilder info= new StringBuilder();
//                   trsctInfoMap.put(trsctList.get(0),trsctInfoList);
                   for(int i=0;i<trsctInfoList.size();i++){
                       List<String> mList = PStringUtils.fastSplit(trsctInfoList.get(i));
                   tems = mList.toArray(new String[mList.size()]);
               info= new StringBuilder();
                     info.append("gene_biotype \"protein_coding\"").append(";").append(tems[8]);
                        if(tems[2].equals("exon")){
                      for(int j=0;j<tems.length-1;j++){
                     bw.write(tems[j]);
                     bw.write("\t");
                     }
                     bw.write(info.toString());
                     bw.newLine();}
                     
                    if(tems[2].equals("CDS")){
                        cdsPosList.add(Integer.parseInt(tems[3]));
                         cdsPosList.add(Integer.parseInt(tems[4])); 
          
                     }
 
                   }
                   
                     strand=tems[6];
                    Collections.sort(cdsPosList);
               if(strand.equals("+")){
for(int k=0;k<cdsPosList.size()/2;k++){ 
       if(cdsPosList.size()==2)continue;
          bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]); bw.write("\t");bw.write("CDS");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(2*k)));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(2*k+1)));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]);bw.write("\t"); bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();
                     }
 bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]); bw.write("\t");bw.write("CDS");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(cdsPosList.size()-2)));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(cdsPosList.size()-1)-3));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]);bw.write("\t"); bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();     
          bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]); bw.write("\t");bw.write("start_codon");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(0)));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(0)+2));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]);bw.write("\t"); bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();
              bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]);bw.write("\t"); bw.write("stop_codon");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(cdsPosList.size()-1)-2));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(cdsPosList.size()-1)));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]);bw.write("\t"); bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();
                    }
                   if(strand.equals("-")){
         for(int k=1;k<cdsPosList.size()/2;k++){
                 
          bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]); bw.write("\t");bw.write("CDS");bw.write("\t");
          bw.write(String.valueOf(cdsPosList.get(2*k)));
          bw.write("\t");bw.write(String.valueOf(cdsPosList.get(2*k+1)));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]);bw.write("\t"); bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();
                     }
 bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]); bw.write("\t");bw.write("CDS");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(0)+3));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(1)));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]);bw.write("\t"); bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();   
          bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]); bw.write("\t");bw.write("start_codon");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(cdsPosList.size()-1)-2));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(cdsPosList.size()-1)));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]); bw.write("\t");bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();
              bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]);bw.write("\t"); bw.write("stop_codon");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(0)));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(0)+2));bw.write("\t");
          bw.write(tems[5]); bw.write("\t");bw.write(tems[6]); bw.write("\t");bw.write(tems[7]);bw.write("\t");bw.write(info.toString());bw.newLine();
                    } 
cdsPosList=new ArrayList<>();
                   tmp=new ArrayList<>();
                   trsctSet=new HashSet<>();
                   trsctInfoList=new ArrayList<>();
                   trsctList=new ArrayList<>();
                    trsctSet.add(trsctID);
                   trsctList.add(trsctID);
                   tmp.add(temp);
            }
            }
            
            trsctInfoList=tmp;
                  StringBuilder info= new StringBuilder();
               for(int i=0;i<trsctInfoList.size();i++){
                List<String> mList = PStringUtils.fastSplit(trsctInfoList.get(i));
                   tems = mList.toArray(new String[mList.size()]);
             info= new StringBuilder();
                     info.append("gene_biotype \"protein_coding\"").append(";").append(tems[8]);
                     if(tems[2].equals("exon")){
                     for(int j=0;j<tems.length-1;j++){
                     bw.write(tems[j]);
                     bw.write("\t");
                     }
                     bw.write(info.toString());
                     bw.newLine();}
                     
                    if(tems[2].equals("CDS")){
                        cdsPosList.add(Integer.parseInt(tems[3]));
                         cdsPosList.add(Integer.parseInt(tems[4]));   
                     }
                 
                   }
                   strand=tems[6];
                    Collections.sort(cdsPosList);
                  if(strand.equals("+")){
for(int k=0;k<-1+cdsPosList.size()/2;k++){    
    if(cdsPosList.size()==2)continue;
          bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]); bw.write("\t");bw.write("CDS");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(2*k)));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(2*k+1)));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]);bw.write("\t"); bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();
                     }
 bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]); bw.write("\t");bw.write("CDS");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(cdsPosList.size()-2)));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(cdsPosList.size()-1)-3));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]);bw.write("\t"); bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();     
          bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]); bw.write("\t");bw.write("start_codon");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(0)));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(0)+2));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]);bw.write("\t"); bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();
              bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]);bw.write("\t"); bw.write("stop_codon");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(cdsPosList.size()-1)-2));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(cdsPosList.size()-1)));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]);bw.write("\t"); bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();
                    }
                   if(strand.equals("-")){
         for(int k=1;k<cdsPosList.size()/2;k++){
                 
          bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]); bw.write("\t");bw.write("CDS");bw.write("\t");
          bw.write(String.valueOf(cdsPosList.get(2*k)));
          bw.write("\t");bw.write(String.valueOf(cdsPosList.get(2*k+1)));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]);bw.write("\t"); bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();
                     }
 bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]); bw.write("\t");bw.write("CDS");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(0)+3));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(1)));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]);bw.write("\t"); bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();   
          bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]); bw.write("\t");bw.write("start_codon");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(cdsPosList.size()-1)-2));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(cdsPosList.size()-1)));bw.write("\t");
          bw.write(tems[5]);bw.write("\t"); bw.write(tems[6]); bw.write("\t");bw.write(tems[7]);bw.write("\t"); bw.write(info.toString());bw.newLine();
              bw.write(tems[0]);bw.write("\t"); bw.write(tems[1]);bw.write("\t"); bw.write("stop_codon");bw.write("\t");bw.write(String.valueOf(cdsPosList.get(0)));bw.write("\t");bw.write(String.valueOf(cdsPosList.get(0)+2));bw.write("\t");
          bw.write(tems[5]); bw.write("\t");bw.write(tems[6]); bw.write("\t");bw.write(tems[7]);bw.write("\t");bw.write(info.toString());bw.newLine();
                    } 
               br.close();
               bw.flush();
               bw.close();

                } catch (Exception e) {
            System.out.println("Error in convert IWGSC GTF to sift4g GTF! please try again! Come on !");
            e.printStackTrace();
        } 
            System.out.println("Congratulations! IWGSC GTF converts to sift4g GTF sucessfully!");
    }
    
    public static void main (String[]args){
    
    new FormatConvertTools();
    }
}
