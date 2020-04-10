/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import static lipengKang.analysis.RandomArray.randomLongCommon;
import utils.IOUtils;
import utils.PStringUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author kanglipeng
 */
public class sample {

    String temp = null;
    String tem[] = null;

    public sample() {
        this.sampleGERP();
    }

    public void sampleGERP() {
        String x[] = {
            "gerp.bed",
            "phyloP_0.bed",
            "phyloP_1.bed",
            "phastCon.bed",
            "CE.bed",
            
            "GSM3449722_macs_CS_seedlings_H3K4me3_ChIP_seq_peaks.bed",
            "GSM3449724_macs_CS_seedlings_H3K4me1_ChIP_seq_rep1_peaks.bed",
            "GSM3449725_macs_CS_seedlings_H3K4me1_ChIP_seq_rep2_peaks.bed",
            "GSM3449726_macs_CS_seedlings_H3K27ac_ChIP_seq_rep1_peaks.bed",
            "GSM3449727_macs_CS_seedlings_H3K27ac_ChIP_seq_rep2_peaks.bed",
           
            "GSM3449728_macs_CS_seedlings_H3K9me2_ChIP_seq_rep1_peaks.bed",
            "GSM3449729_macs_CS_seedlings_H3K9ac_ChIP_seq_rep1_peaks.bed",
            "GSM3449730_macs_CS_seedlings_H3K9ac_ChIP_seq_rep2_peaks.bed",
            "GSM3449731_macs_CS_seedlings_H3K36me3_ChIP_seq_rep1_peaks.bed",
            "GSM3449732_macs_CS_seedlings_H3K36me3_ChIP_seq_rep2_peaks.bed",
            
            "GSM3564341_macs_CS_seedlings_H3K9me2_ChIP_seq_rep2_peaks.bed",
            "GSM3564342_macs_CS_seedlings_DNaseI_seq_rep1_peaks.bed",
            "GSM3564343_macs_CS_seedlings_DNaseI_seq_rep2_peaks.bed",
            "GSM3901000_CS_seedlings_lncRNA_seq_rep1.rpm.bed",
            "GSM3901002_CS_seedlings_lncRNA_seq_rep2.rpm.bed",
            
            "GSM4043024_macs_CS_seedlings_H3K27me3_ChIP_seq_peaks.bed",
            "GSM3449734_CS_seedlings_Bisulfite_seq.CG.bed",
            "GSM3449734_CS_seedlings_Bisulfite_seq.CHG.bed",
            "GSM3449734_CS_seedlings_Bisulfite_seq.CHH.bed",    
        };

        // for(String j:y){
/////    System.out.println("zcat "+j+".gz | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5}' >"+j);
////    }
//System.out.print("unionBedGraphs -k ");
//for (String k:x){
////     System.out.println( "grep -v '#bedGraph' "+k+" | awk '{print $4}' >x.bed");
////     System.out.println( "grep -v '#bedGraph' "+k+" | awk '{print $1\"\\t\"$2\"\\t\"$3}' >y.bed");
////  System.out.println("awk '{if($1 ~ /^[0-9]+$/) print $1;else printf \"%2.3f\\n\",$1;}' x.bed >xx.bed");
//// System.out.println("paste -d \"\\t\" y.bed xx.bed >"+k+".bed");
//
//System.out.print(k+" ");
//
//}System.out.print(">genomic_features.ins");
//System.out.println("java -jar PlantGenetics.jar -m /data1/home/lipeng/database/epigenomic/"+k+" -o /data1/home/lipeng/database/epigenomic/"+k+".bed");
//System.out.println("mv "+k+".bed "+k);
//System.out.print(k+" ");
//    }
//  System.out.print(">genomic_features.bedGraph");  
//  int chrA[]={471304005,462376173,454103970,452555092,453230519,452440856,450046986};
//        int chrB[]={438720154,453218924,448155269,451014251,451372872,452077197,453822637};
//        int chrD[]={452179604,462216879,476235359,451004620,451901030,450509124,453812268};
        String Genome[] = {"A", "B", "D"};
       // String block[] = {"0d", "3UTR", "5UTR", "cns", "intergenic", "intron"};
          String block[] = {"3UTR", "5UTR", "cds", "geneFlank2kb","intron"};
        int temp[] = {19231931, 4556208, 2440904, 4071004, 128887516, 32740874,
            19045278, 4466587, 2438867, 4011128, 125193507, 32194716,
            19625114, 4705418, 2569894, 4213555, 130350911, 32816806};

//  String querySpecies[]={"ancom","ortho","ercur","ertef","sobic","zemay","eccru", "seita", "pahal","pavir","phedu","brdis","hovul","hospo","aetau","traes"+j,"trura","orbra","leper","orpun","ormer","orlon","orglu","orbar","orgla","orind","orniv","orruf","orjap"};    
        try {

            for (int i = 1; i <= 1; i++) {
                String outDir = "/Users/kanglipeng/Desktop/1A/" + i + ".sh";
                BufferedWriter bw = IOUtils.getTextWriter(outDir);
                int mark = 0;
                   for (String j : Genome) {
               for (int m = 1; m <= 7; m++) { 
             // System.out.println("grep  -f "+m+j+".longestID.txt "+m+j+".bed >"+m+j+".long.bed"); 
                System.out.println("grep  -f "+m+j+".longestID.txt "+m+j+".bed >"+m+j+".long.bed");  
//                       System.out.println("zcat "+m+j+".INSIGHTs.gz |sed 's/chr"+m+"/chr"+m+j+"/g' >"+m+j+".INSIGHT "); 
//                       System.out.println("/data1/home/lipeng/software/LINSIGHT/LINSIGHT-prep -i "+m+j+".INSIGHT -o "+j+"/ "); 
  //for (String k : x) {
                        // System.out.println("/data1/home/lipeng/software/est-sfs-release-2.03/est-sfs /data1/home/lipeng/software/est-sfs-release-2.03/config-JC.txt /data2/lipeng/asAlle/est/est-usfs/"+j+"."+k+".input /data1/home/lipeng/software/est-sfs-release-2.03/seedfile.txt /data2/lipeng/asAlle/est/est-usfs/"+j+"."+k+".usfs /data2/lipeng/asAlle/est/"+j+"."+k+".probs ");
                        //  System.out.println("grep '"+j+"' ../../../ncDelSNP/"+k+".bed | awk '{print substr($1,1,4)\"\\t\"$2\"\\t\"$3}' >"+j+"."+k+".bed");
                        //  System.out.println("bedtools intersect -a "+j+"SH.flt.est -b "+j+"."+k+".bed >"+j+"."+k+".est");
//     
 //System.out.println("bedtools subtract -a  "+m+j+"."+k+" -b cds.bed >"+m+j+".nc."+k);
// System.out.print(m+j+"."+k+" ");

 //System.out.println("bedtools intersect -a "+j+".vcf -b "+j+"."+k+".bed >"+j+"."+k+".vcf");
 //System.out.println("awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' "+j+"."+k+".bed >>xx");

//##########################dfe-est#######################
//                        bw.write("shuf -n " + temp[mark] / 2 + " /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + ".est|sortBed -i >/data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".temp_0");
//                        bw.newLine();
//                        bw.write("cut -f4-6 /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".temp_0 >/data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".input");
//                        bw.newLine();
//                        bw.write("/data1/home/lipeng/software/est-sfs-release-2.03/est-sfs /data1/home/lipeng/software/est-sfs-release-2.03/config-JC.txt /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".input /data1/home/lipeng/software/est-sfs-release-2.03/seedfile.txt /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".usfs /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k +"." + i + ".probs ");
//                        bw.newLine();
//                        bw.write("cat /data2/lipeng/asAlle/est/est-usfs/head  /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".usfs /data2/lipeng/asAlle/est/est-usfs/" + j + ".4d.usfs|sed 's/,/ /g' >" + j + "." + k + "." + i + ".usfs");
//                        bw.newLine();
//                        bw.write("sed 's/A.3UTR/" + j + "." + k + "." + i + "/g' example-config-file-for-est_dfe-site_class-0.txt >" + j + "." + k + "." + i + "-0.mat");
//                        bw.newLine();
//                        bw.write("./est_dfe -c " + j + "." + k + "." + i + "-0.mat");
//                        bw.newLine();
//                        bw.write("sed 's/A.3UTR/" + j + "." + k + "." + i + "/g' example-config-file-for-est_dfe-site_class-1.txt >" + j + "." + k + "." + i + "-1.mat");
//                        bw.newLine();
//                        bw.write("bedtools subtract -a /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".temp_0 -b /data2/lipeng/asAlle/est/est-usfs/" + j + "SH.rare.est >/data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".temp_1");
//                        bw.newLine();
//                        bw.write(" grep -e '100,' -e ',100'  /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".temp_1 |cut -f5 |awk -F ',' '{if($1==1)print \"1\" ;else if($2==1)print \"2\";else if($3==1)print \"3\";else if($4==1)print \"4\" }'  >" + j + "." + k + "." + i + ".temp_2");
//                        bw.newLine();
//                        bw.write(" grep -e '100,' -e ',100'  /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".temp_1| cut -f4| awk -F ',' '{if($1==100)print \"1\" ;else if($2==100)print \"2\";else if($3==100)print \"3\";else if($4==100)print \"4\" }'  >" + j + "." + k + "." + i + ".temp_1");
//                        bw.newLine();
//                        bw.write(" grep -e '100,' -e ',100'  /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".temp_1|cut -f6|awk -F ',' '{if($1==1)print \"1\" ;else if($2==1)print \"2\";else if($3==1)print \"3\";else if($4==1)print \"4\" }'  >" + j + "." + k + "." + i + ".temp_3");
//                        bw.newLine();
//                        bw.write("paste " + j + "." + k + "." + i + ".temp_1 " + j + "." + k + "." + i + ".temp_2 " + j + "." + k + "." + i + ".temp_3" + "|awk '{if($1!=$2&&$1!=$3)print $0}'|wc -l >" + j + "." + k + "." + i + ".temp_0");
//                        bw.newLine();
//                        bw.write("wc -l /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k +"." + i + ".input|awk '{print $1}'>>" + j + "." + k + "." + i + ".temp_0");
//                        bw.newLine();
//                        bw.write("awk 'NR==1{a=$0}NR==2{print a/$0}' " + j + "." + k + "." + i + ".temp_0 |awk '{print \"p_additional\"\"\\t\"$0\"\\n\"\"s_additional\"\"\\t\"\"0.1\"}' >>" + j + "." + k + "." + i + "-1.mat");
//                        bw.newLine();
//                        bw.write("./est_dfe -c " + j + "." + k + "." + i + "-1.mat");
//                        bw.newLine();
//                        bw.write("./prop_muts_in_s_ranges -c " + j + "." + k + "." + i + "_sel/est_dfe.out -o prop_muts/" + j + "." + k + "." + i + ".prop_muts");
//                        bw.newLine();
//                        bw.write("sed 's/A.3UTR/" + j + "." + k + "." + i + "/g' example-config-file-for-est_alpha_omega.txt |sed  's/ est_alpha_omega.out/alpha_omega\\/" + j + "." + k + "." + i + ".alpha_omega.out/g' >" + j + "." + k + "." + i + "-2.mat");
//                        bw.newLine();
//                        bw.write("awk ' BEGIN { ORS=\" \" } { print }' " + j + "." + k + "." + i + ".temp_0 |awk 'BEGIN {FS=\" \"; OFS=\" \"} {print $2, $1}'|awk '{print \"1 \"$1\" \"$2}' >" + j + "." + k + "." + i + ".divergent.txt");
//                        bw.newLine();
//                        bw.write("awk ' BEGIN { ORS=\" \" } { print }' " + j + ".4d.temp |awk 'BEGIN {FS=\" \"; OFS=\" \"} {print $2, $1}'|awk '{print \"0 \"$1\" \"$2}' >>" + j + "." + k + "." + i + ".divergent.txt");
//                        bw.newLine();
//                        bw.write("./est_alpha_omega -c " + j + "." + k + "." + i + "-2.mat");
//                        bw.newLine();
//                        bw.write("rm -rf " + j + "." + k + "." + i + "_neut");
//                        bw.newLine();
//                         bw.write("rm  -rf " + j + "." + k + "." + i + "_sel");
//                        bw.newLine();
//                        bw.write("rm " + j + "." + k + "." + i + ".temp_3");
//                        bw.newLine();
//                        bw.write("rm " + j + "." + k + "." + i + ".temp_1");
//                        bw.newLine();
//                        bw.write("rm " + j + "." + k + "." + i + ".temp_0");
//                        bw.newLine();   
//                         bw.write("rm " + j + "." + k + "." + i + ".temp_2");
//                        bw.newLine();  
//                         bw.write("rm " + j + "." + k + "." + i + ".divergent.txt");
//                        bw.newLine();
//                        bw.write("rm /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".temp_1");
//                        bw.newLine();
//                        bw.write("rm /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".input");
//                        bw.newLine();
//                           bw.write("rm /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".probs");
//                        bw.newLine();
//                         bw.write("rm /data2/lipeng/asAlle/est/est-usfs/" + j + "." + k + "." + i + ".temp_0");
//                        bw.newLine();
//                        bw.write("rm "+ j + "." + k + "." + i + "-1.mat");
//                        bw.newLine();
//                         bw.write("rm "+ j + "." + k + "." + i + "-2.mat");
//                        bw.newLine();
//                         bw.write("rm "+ j + "." + k + "." + i + "-0.mat");
//                        bw.newLine();
//                        bw.write("rm /data2/lipeng/asAlle/est/est-usfs/"+ j + "." + k + "." + i + ".usfs");
//                        bw.newLine();
//                        bw.write("rm "+ j + "." + k + "." + i + ".usfs");
//                        bw.newLine();
//                        mark++;

                        
                        
                        
                        
                        
     //  System.out.println("for i in {0..19};do cat software/dfe-alpha-release-2.16-${i}/prop_muts/"+j+"."+k+"* >>"+j+"."+k+".prop_muts & done");
                        //  System.out.println("bedtools intersect -b "+i+j+".Dnase_seq -a "+i+j+".Nucleosome >"+i+j+".NDR");
//System.out.println("java -jar PlantGenetics.jar -m /data2/lipeng/asAlle/temp/test/"+i+j+".maf -o /data2/lipeng/asAlle/temp/test/"+i+j+".fcoal.est -g traes"+j);
//System.out.println("java -jar dfe_outGroup.jar -m /data2/lipeng/asAlle/temp/test/"+i+j+".maf -o /data2/lipeng/asAlle/temp/test/"+i+j+".outGrop.est -g traes"+j);
//System.out.println("bedtools intersect -a /data2/lipeng/asAlle/temp/test/"+i+j+".outGrop.est -b /data2/lipeng/asAlle/temp/test/"+i+j+".focal.est >/data2/lipeng/asAlle/temp/test/temp1");
//System.out.println("bedtools intersect -b /data2/lipeng/asAlle/temp/test/"+i+j+".outGrop.est -a /data2/lipeng/asAlle/temp/test/"+i+j+".focal.est >/data2/lipeng/asAlle/temp/test/temp2");
//System.out.println("paste /data2/lipeng/asAlle/temp/test/temp2  /data2/lipeng/asAlle/temp/test/temp1 | cut -f1-4,8-9 >/data2/lipeng/asAlle/temp/test/"+i+j+".est");
//System.out.print("cat "+j+"."+k+".out ");
//System.out.println("sed -i 's/,/ /g' "+j+"."+k+".est");
//System.out.println("sed -i 's/"+j+"."+k+"_results_file/est_alpha_omega_results_file/g' "+j+"."+k+"-2.mat");
//System.out.println("sed -i 's/"+i+j+"_2/"+i+j+"/g' "+i+j+".cons.bed");
//##usfs generator
//System.out.println("./ml-est-sfs-two-outgroup "+j+"."+k+".est "+j+"."+k+".usfs");
//
//##dfe config mat generator
//System.out.println("sed 's/A.3UTR/"+j+"."+k+"/g' example-config-file-for-est_dfe-site_class-0.txt >"+j+"."+k+"-0.mat");
//System.out.println("sed 's/A.3UTR/"+j+"."+k+"/g' example-config-file-for-est_dfe-site_class-1.txt >"+j+"."+k+"-1.mat");
//System.out.println("sed 's/A.3UTR/"+j+"."+k+"/g' example-config-file-for-est_alpha_omega.txt  >"+j+"."+k+"-2.mat");
//##run dfe
//System.out.println("./est_dfe -c "+j+"."+k+"-0.mat");
//System.out.println("./est_dfe -c "+j+"."+k+"-1.mat");
//System.out.println("./est_alpha_omega -c "+j+"."+k+"-2.mat");
//System.out.println("./prop_muts_in_s_ranges -c "+j+"."+k+"_sel/est_dfe.out -o "+j+"."+k+".prop_file");
//################dfe-config-1 P_additional and s_additional adder
//System.out.println("grep '100' "+j+"."+k+".est|cut -f2|awk '{if($1==1)print \"1\" ;else if($2==1)print \"2\";else if($3==1)print \"3\";else if($4==1)print \"4\" }'  >x");
//System.out.println("grep '100' "+j+"."+k+".est|awk '{if($1==100)print \"1\" ;else if($2==100)print \"2\";else if($3==100)print \"3\";else if($4==100)print \"4\" }'  >xx");
//System.out.println("paste x xx |awk '{if($1!=$2)print $0}'|wc -l >"+j+"."+k+".class1");
//System.out.println("grep '100' "+j+"."+k+".est|wc -l >>"+j+"."+k+".class1");
//System.out.println("awk 'NR==1{a=$0}NR==2{print a/$0}' "+j+"."+k+".class1 |awk '{print \"p_additional\"\"\\t\"$0\"\\n\"\"s_additional\"\"\\t\"\"0.1\"}' >>../dfe-alpha-release-2.16/"+j+"."+k+"-1.mat");
////
//#################dfe-diverget.txt must used after dfe-config-1
//System.out.println("wc -l "+j+"."+k+".est|awk '{print \"1 \"$1}' >x");
//System.out.println("head -n 1 "+j+"."+k+".class1|paste -d \" \" x - >>../dfe-alpha-release-2.16/"+j+"."+k+".divergent.txt");
//System.out.println("wc -l "+j+".4d_sites.est|awk '{print \"0 \"$1}' >x");
//System.out.println("head -n 1 "+j+".4d_sites.class1|paste -d \" \" x - >>../dfe-alpha-release-2.16/"+j+"."+k+".divergent.txt");
//################usfs 
//System.out.println("awk '{print \"1\"\"\\n\"\"100\"}' x >../dfe-alpha-release-2.16/"+j+"."+k+".usfs");
//System.out.println("cat "+j+"."+k+".usfs "+j+".4d_sites.usfs >>../dfe-alpha-release-2.16/"+j+"."+k+".usfs" );
//System.out.println("sed -i 's/,/ /g' "+j+"."+k+".usfs");
// System.out.println("bedtools intersect -a "+i+j+".bed -b ../phastCon/"+j+".phastCon.sample >"+i+j+".gerp.sample");
//  System.out.println("paste "+i+j+"._02.est "+i+j+"._01.est |cut -f1-4,8-9 >"+i+j+"._01.input");
                        //System.out.println("rm "+i+j+"_2.bg");
                        // System.out.println("rm "+i+j+"_1.bg");
// System.out.println("tail -n +9 "+i+j+".anc |awk '{print $3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' >x");
// System.out.println("paste ../../"+i+j+"._1txt x >../../"+i+j+"._2txt");
//System.out.println("bedtools subtract -a "+i+j+".input -b "+i+j+"._01.input >"+i+j+".temp");
//System.out.println("cut -f 4-6 "+i+j+".temp |sed 's/,/ /g'| sed 's/20/100/g' >"+i+j+".temp1");
//System.out.println("cut -f1-3 "+i+j+".temp >"+i+j+".temp2");
//System.out.println("paste "+i+j+".temp2 "+i+j+".temp1 >"+i+j+".temp3");
//System.out.println("sed 's/,/ /g' "+i+j+"._01.input >" +i+j+".temp2");
//System.out.println("cat "+i+j+".temp3 "+i+j+".temp2 >"+i+j+".temp1");
//System.out.println("sortBed -i "+i+j+".temp1 >"+i+j+".est");
//System.out.println("rm "+i+j+".temp1");System.out.println("rm "+i+j+".temp2");System.out.println("rm "+i+j+".temp");System.out.println("rm "+i+j+".temp3");
//  System.out.println("cut -f1 "+j+"."+k+".usfs.input|awk -F' ' '{printf(\"%.f \",$1/5);printf(\"%.f \",$2/5);printf(\"%.f \",$3/5);printf(\"%.f\\n\",$4/5)}' >x");
//System.out.println("cut -f1 "+i+j+".input|awk -F',' '{printf(\"%.f,\",$1/5);printf(\"%.f,\",$2/5);printf(\"%.f,\",$3/5);printf(\"%.f\\n\",$4/5)}' >x");
//System.out.println("cut -f2-3 "+i+j+".input >xx");
//System.out.println("paste x xx >" +i+j+".usfs.input");
//System.out.println("/data1/home/lipeng/software/est-sfs-release-2.03/est-sfs /data1/home/lipeng/software/est-sfs-release-2.03/config-rate6_1.txt "+i+j+".usfs.input /data1/home/lipeng/software/est-sfs-release-2.03/seedfile.txt "+i+j+".usfs "+i+j+".anc");
//System.out.println("cut -f2-3 "+j+"."+k+".usfs.input >xx");
//System.out.println("paste x xx >" +j+"."+k+".usfs.input");
//System.out.println("./ml-est-sfs-two-outgroup " +j+"."+k+".usfs.input "+j+"."+k+".usfs.output");
//System.out.println("awk '{print substr($1,1,4)\"\\t\"$2\"\\t\"$3}'  /data2/lipeng/ncDelSNP/"+j+"."+k+".bed >/data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".bed");
//System.out.println("mafsInRegion /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".bed /data2/lipeng/asAlle/est/temp/A/"+ j+"."+k+".snp.maf"+" /data2/lipeng/asAlle/est/temp/A/"+j+"SH.flt.maf");
//System.out.println("mafsInRegion /data2/lipeng/msa/multiz"+j+"/4d_"+j+"_sites.bed /data2/lipeng/asAlle/est/temp/A/"+ j+"."+k+".snp.maf"+" /data2/lipeng/asAlle/est/temp/A/"+j+"SH.flt.maf");
//System.out.println("grep -v \"null\" /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp.maf |awk '{if(substr($7,1,1)!=\"-\")print $0}' >/data2/lipeng/asAlle/est/temp/A/"+j+"."+k+"_00.bed");
//System.out.println("mafFilter -minRow=3 /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+"_00.bed >/data2/lipeng/asAlle/est/temp/A/"+ j+"."+k+"_01.bed");
//System.out.println("java -jar PlantGenetics.jar -m /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+"_01.bed "+"-o /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp.txt"+" -g traes"+j);
//System.out.println("bedtools intersect -b /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp.txt"+" -a /data2/lipeng/asAlle/est/temp/A/"+j+".focal.est >/data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp_0.txt");
// System.out.println("bedtools intersect -a /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp.txt"+" -b /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp_0.txt >/data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp_1.txt");
//  System.out.println("sortBed -i /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp_0.txt >/data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp_2.txt");
//    System.out.println("sortBed -i /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp_1.txt >/data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp_3.txt");
//  System.out.println("paste /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp_2.txt /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp_3.txt >/data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp_4.txt"); 
//  System.out.println("cut -f 1-4,8,9 /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".snp_4.txt  >/data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".usfs.input"); 
//  System.out.println("rm /data2/lipeng/asAlle/est/temp/A/*.snp.maf");
//    System.out.println("rm /data2/lipeng/asAlle/est/temp/A/*.txt");
//     System.out.println("rm /data2/lipeng/asAlle/est/temp/A/*.bed");
//     System.out.println("cut -f4-6 /data2/lipeng/asAlle/est/temp/A/"+j+"."+k+".usfs.input |sed 's/,/ /g' > /data1/home/lipeng/software/estimate-unfolded-sfs-release-1.02/"+j+"."+k+".usfs.input");
//     System.out.println("/data1/home/lipeng/software/estimate-unfolded-sfs-release-1.02/ml-est-sfs-two-outgroup /data1/home/lipeng/software/estimate-unfolded-sfs-release-1.02/"+j+"."+k+".usfs.input /data1/home/lipeng/software/estimate-unfolded-sfs-release-1.02/"+j+"."+k+".usfs.output");
                        // String outDir = "/Users/kanglipeng/Desktop/1A/" + +i + j + ".sh";
                        // BufferedWriter bw = IOUtils.getTextWriter(outDir);
//                   
//                    System.out.println("paste -a "+j+"x -b "+j+"s >"+j+"m");
//                     System.out.println("cut -f1-4,8-10 "+j+"mm >traes"+j+".hovul.brdis.phedu.est");
//                     System.out.println("rm "+j+"x");
//                     System.out.println("rm "+j+"m");
//                     System.out.println("rm "+j+"s");
//                     System.out.println("rm "+j);
//System.out.println("nohup convert2bed -i wig <"+i+j+".phyloP.wig> "+i+j+".phyloP.bed &");
                        //System.out.println("awk '{printf(\"%s\\t %s\\t %s\\t %.3f\\n\", $1,$2,$3,$5)}' "+i+j+".phyloP.bed >"+i+j+".phyloP.beds");//bw.newLine();
                        //System.out.println("mv "+i+j+".phyloP.beds "+i+j+".phyloP.bed");
                        //System.out.println("nohup gunzip "+i+j+".phastCon.wig.gz &");//}
                        //     StringBuilder all = new StringBuilder();
                        //   for (String m : x) {
                        //        all.append(" " + i + j + "." + m);
                        //   }
                        // bw.write("bedtools unionbedg -i" + all+" >"+i+j+".features.bed");
//  bw.write("head -n 800000 temp/"+i+j+".bed >"+i+j+"_part1.bed");bw.newLine();
//    bw.write("tail -n +800001 temp/"+i+j+".bed >"+i+j+"_part2.bed");bw.newLine();
// bw.write("mafsInRegion "+i+j+"_part1.bed "+i+j+"_part1.maf temp/"+i+j+".sort.maf");bw.newLine();
// bw.write("mafsInRegion "+i+j+"_part2.bed "+i+j+"_part2.maf temp/"+i+j+".sort.maf");bw.newLine();
//  bw.write("phastCons --target-coverage 0.25 --expected-length 12 --rho 0.4 --most-conserved "+i+j+"_part1.CE.bed "+i+j+"_part1.maf ../phylop/"+j+".mod > "+i+j+"_part1.phastCon.wig 2>x");bw.newLine();
//bw.write("phastCons --target-coverage 0.25 --expected-length 12 --rho 0.4 --most-conserved "+i+j+"_part2.CE.bed "+i+j+"_part2.maf ../phylop/"+j+".mod > "+i+j+"_part2.phastCon.wig 2>x");bw.newLine();
//    bw.write("cat "+i+j+"_part1.CE.bed "+i+j+"_part2.CE.bed >"+i+j+".CE.bed");bw.newLine();
// bw.write("rm "+i+j+"_part1.CE.bed");bw.newLine();
// bw.write("rm "+i+j+"_part2.CE.bed");bw.newLine();
//    bw.write("cat "+i+j+"_part1.phastCon.wig "+i+j+"_part2.phastCon.wig >"+i+j+".phastCon.wig");bw.newLine();
//bw.write("rm "+i+j+"_part1.phastCon.wig");bw.newLine();
//    bw.write("rm "+i+j+"_part2.phastCon.wig");bw.newLine();
//bw.write("sed 's/_part1//g' "+i+j+".CE.bed |sed 's/_part2//g' >"+i+j+".CE.beds");bw.newLine();
//bw.write("mv "+i+j+".CE.beds "+i+j+".CE.bed");bw.newLine();
//bw.write("sed 's/_part1//g' "+i+j+".phastCon.wig |sed 's/_part2//g' > "+i+j+".phastCon.wigs");bw.newLine();
//bw.write("mv "+i+j+".phastCon.wigs "+i+j+".phastCon.wig");
//System.out.println("mafsInRegion x.bed "+k+j+".maf "+j+".simple.maf");
//System.out.println("");
//System.out.println("prequel -k -r /data1/home/lipeng/database/genome/traes"+j+"/"+k+j+".fa -k MAF "+k+j+".maf "+j+".mod "+k+j);
//System.out.println("mkdir "+i+j );
//System.out.println("grep 'chr"+i+"'  ../../asGenome/"+j+"/split/"+j+".txt >"+i+j+"/x.bed" );
//System.out.println("cd "+i+j );
//for(int k=1;k<=80;k++){
//System.out.println("sed -n '"+k+"p' x.bed >x" );
////System.out.println("mafsInRegion x "+k+".maf ../"+j+".maf");
//System.out.println("bedtools intersect -a x -b ../"+j+".pns.bed >"+k+".bed");
//System.out.println("rm x");}
//System.out.println("cd ..");
//System.out.println("awk '{print $1\"\\t\"0\"\\t\"$2}' ../../msa/gerp/"+i+j+".fai >"+i+j+".bed");
//System.out.println("mafsInRegion "+i+j+".bed "+i+j+".maf 27way."+j+".maf &");
//System.out.println("java -jar gerpReformat.jar -f /data2/lipeng/asAlle/GERP/"+i+j+".bed -m /data2/lipeng/asAlle/GERP/"+i+j+".maf.rates -o /data2/lipeng/asAlle/GERP/"+i+j+".wig -g "+j);
//System.out.println("grep 'chr"+i+j+"' /data1/home/lipeng/database/genome/traes/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta.fai|awk '{print substr($1,1,4)\"\\t\"$2}' >"+i+j+".fai");
//System.out.println("multiBigwigSummary bins -b "+i+j+".gerp.bw -o results.npz -bs 1000000 --outRawCounts "+i+j+".tab");
//System.out.println("zcat /data1/publicData/wheat/genotype/VMapII/VMap2.1/chr0"+((i-1)*6+5)+"_vmap2.1.vcf.gz |tail -n +22|awk '{split($8,a,\";\");print \"chr"+i+"\"\"\\t\"$2-1\"\\t\"$2\"\\t\"$4\"\\t\"$5\"\\t\"substr(a[7],5,length(a[7]))\"\\t\"substr(a[8],9,length(a[8]))\"\\t\"substr(a[9],8,length(a[9]))}' >/data2/lipeng/work/"+i+j+"_part1.vcf ");
//  System.out.println("zcat /data1/publicData/wheat/genotype/VMapII/VMap2.1/chr0"+((i-1)*6+6)+"_vmap2.1.vcf.gz |tail -n +22|awk '{split($8,a,\";\");print \"chr"+i+"\"\"\\t\"$2+"+(chrD[i-1]-1)+"\"\\t\"$2+"+chrD[i-1]+"\"\\t\"$4\"\\t\"$5\"\\t\"substr(a[7],5,length(a[7]))\"\\t\"substr(a[8],9,length(a[8]))\"\\t\"substr(a[9],8,length(a[9]))}' >/data2/lipeng/work/"+i+j+"_part2.vcf ");
                        // System.out.println("cat "+i+j+"_part1.vcf "+i+j+"_part2.vcf >"+i+j+".vcf");
//bw.write("awk '{print $1\"\\t\"0\"\\t\"$2}' /data2/lipeng/msa/gerp/"+i+j+".fai >/data2/lipeng/asAlle/maf/"+i+j+".bed");bw.newLine();
//bw.write("mafsInRegion /data2/lipeng/asAlle/maf/"+i+j+".bed /data2/lipeng/asAlle/maf/"+i+j+".maf /data2/lipeng/asAlle/maf/27way."+j+".maf ");bw.newLine();
//bw.write("mafRanges /data2/lipeng/asAlle/maf/" +i+j+".maf traes"+j+" /data2/lipeng/asAlle/maf/"+i+j+".bed");bw.newLine();
//bw.write("gerpcol -f /data2/lipeng/asAlle/maf/"+i+j+".maf -t /data2/lipeng/asAlle/tree/27way."+j+".tree -e traes"+j+" -j -z");bw.newLine();
//bw.write("java -jar gerpReformat.jar -f /data2/lipeng/asAlle/maf/" +i+j+".bed -m /data2/lipeng/asAlle/maf/"+i+j+" -o /data2/lipeng/asAlle/maf/"+i+j+".wig");bw.newLine();
//bw.write("bgzip /data2/lipeng/asAlle/maf/"+i+j+".wig");bw.newLine();
                           bw.flush();
                        //        bw.close();
//
//System.out.println("gerpcol -f ../../asAlle/maf/"+i+j+".maf -t ../../asAlle/tree/27way."+j+".tree -e traes"+j+" -j -z");
//System.out.println("sed 's/A/"+j+"/g' "+i+j+".wig >"+i+j+".wigs");
//System.out.println("mv "+i+j+".wigs "+i+j+".wig");
                  //  }
                // System.out.println("-g wheat.genome >"+m+j+".features.bed");  
                    }
                }
                bw.close();

            }

        } catch (Exception e) {
            System.out.println("Error in calling sample GERP!");
            e.printStackTrace();

        }
    }

    public static void main(String[] args) {
        new sample();
    }
}
