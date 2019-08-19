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

    String temp = null;
    String[] tem = null;

    public FormatConvertTools() {
       // this.debugGERP();
      //  this.mafToFa();
      this.gerpNorm();
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
          
            int gerpPos=0;
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
 
 gerpPos=gerpPos+tem[6].replaceAll("N", "").replaceAll("n", "").replaceAll("-", "").length();
 if(tem[2].equals("508641235"))break;
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
        String maf = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.maf";
        String out = "/Users/kanglipeng/Desktop/iwgsc_refseqv1.1_genes_2017July06/test.out";
        BufferedReader br = YaoIOUtils.getTextReader(maf);
        BufferedWriter bw = YaoIOUtils.getTextWriter(out);
        String[]species={"wheatA","wheatD","wheatB","barley"};
         HashMap<String,String> blockInfo=new HashMap();
         StringBuilder absentSeq=new StringBuilder();
        ArrayList<String>[]specieSeq=new ArrayList[species.length];
         for(int l=0;l<species.length;l++){specieSeq[l]=new ArrayList<String>();}
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("s")) {
                    List<String> fList = KStringUtils.fastSplit(temp);
                    List<String> fListNew = new ArrayList<>();
                    for (int i = 0; i < fList.size(); i++) {
                        if (fList.get(i) != null && !fList.get(i).equals("")) {
                            fListNew.add(fList.get(i));
                        }
                    }
 tem = fListNew.toArray(new String[fListNew.size()]);
 blockInfo.put(tem[1].split("\\.")[0],tem[6]);
                }
                
                if(temp.equals("")){
                    ArrayList<String>blockSpecies=new ArrayList(blockInfo.keySet());
                    for(int i=0;i<species.length;i++){
               if(Collections.binarySearch(blockSpecies, species[i])<0){
                    for (int j = 0; j < blockInfo.get(species[0]).length(); j++) {
                                                    absentSeq.append('-');
                                                }
               blockInfo.put(species[i], absentSeq.toString());
               absentSeq=new StringBuilder();
               }
                  specieSeq[i].add(blockInfo.get(species[i]));
                    }
                    blockInfo=new HashMap();
                }
            }
            br.close();
            for(int k=0;k<specieSeq.length;k++){
                bw.write(">"+species[k]);
                bw.newLine();
            for(String m:specieSeq[k]){
            bw.write(m);
            bw.newLine();
            }
            bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in convert IWGSC GTF to sift4g GTF! please try again! Come on !");
            e.printStackTrace();
        }
        System.out.println("Congratulations! IWGSC GTF converts to sift4g GTF sucessfully!");
    }

    //gerpNorm: parse GERP.rates to GERP.fa
    //GERP.rates follows align.maf blocks---->GERP.fa follows chromosome blocks
    public void gerpNorm(){
     String maf = "/data1/home/lipeng/database/genome/barley/Hordeum_vulgare.IBSC_v2.43.longestTrans.gff3";
        String out1 = "/data1/home/lipeng/database/genome/barley/CNEfinderGene.gff3";
        String out2 = "/data1/home/lipeng/database/genome/barley/CNEfinderExon.gff3";
        BufferedReader br = YaoIOUtils.getTextReader(maf);
        BufferedWriter bw1 = YaoIOUtils.getTextWriter(out1);
        BufferedWriter bw2 = YaoIOUtils.getTextWriter(out2);
       
       try {
            bw1.write("gene_name");bw1.write("\t");bw1.write("chromosome_name");bw1.write("\t");bw1.write("start_position");bw1.write("\t");bw1.write("end_position");
            bw1.newLine();
            
            while ((temp = br.readLine()) != null) {
                
                  List<String> tList = PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                  
if(tem[2].equals("gene")){
    bw1.write(tem[8].split(";")[0].split(":")[1]);bw1.write("\t");bw1.write(tem[0]);bw1.write("\t");bw1.write(tem[3]);bw1.write("\t");bw1.write(tem[4]);
    bw1.newLine();
}
  if(tem[2].equals("exon")){
    bw2.write(tem[0]);bw2.write("\t");bw2.write(tem[3]);bw2.write("\t");bw2.write(tem[4]);;bw2.write("\t");bw2.write(tem[8].split(";")[0].split(":")[1]);
    bw2.write("\t");bw2.write(tem[7]);bw2.write("\t");bw2.write(tem[6]);bw2.newLine();
  }          
            
            }br.close();
            bw1.flush();bw1.close();
            bw2.flush();bw2.close();

       } catch (Exception e) {
            System.out.println("Error in convert IWGSC GTF to sift4g GTF! please try again! Come on !");
            e.printStackTrace();
        }
    }
    
    
    public static void main(String[] args) {
        new FormatConvertTools();
    }
}
