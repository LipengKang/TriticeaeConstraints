/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import java.io.BufferedWriter;
import utils.IOUtils;

/**
 *
 * @author kanglipeng
 */
public class CommandPrinter {
//    (((traesD:0.00148569993756236178,aetau:0.00137597763276625803):0.02738384320412087097,trura:0.03299826365421796442):0.02391964411027163590,(hovul:0.00321937415321173637,hospo:0.00299195432052477709):0.05768654924925382954)
//    ((traesB:0.03148792133972425078,(trura:0.03456613721152513446,aetau:0.02701086621543083433):0.00652464237566866608):0.02265886047553810484,(hovul:0.00316931588925377199,hospo:0.00292807003093157421):0.05848959447956827984)
//    ((hovul:0.00315531931731526452,hospo:0.00291301039210765597):0.05739904630277179592,(aetau:0.02889020081603003998,(traesA:0.00385802966340279994,trura:0.00434378444463454777):0.02797814516990456898):0.02394215872622542166)

    public static void main(String[] args) {
//        try {
            //create blocks for all chr
            int chr=7;
            String genome="D";
//            for (int j =1 ; j <= 9; j++) {
//                System.out.println("zcat /data1/publicData/wheat/genotype/VMapII/VMap2.1/chr00"+j+"_vmap2.1.vcf.gz |tail -n +22|awk '{split($8,a,\";\");print \"chr"+j+"\"\"\\t\"$2-1\"\\t\"$2\"\\t\"$4\"\\t\"$5\"\\t\"substr(a[7],5,length(a[7]))}' >chr"+j+"_vmap2.1.vcf.bed");
//                
////                    bw.newLine();j
//                      }
             for (int j =1 ; j <= 7; j++) {
                System.out.println("awk 'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{if($1==\""+j+"A\")$1=\"chr\"substr($1,1,1);print $0 }' vmap2.1.vcf.bed >"+j+"A_vmap2.1.vcf.bed");
                      }
             

              // System.out.println("cat *_vmap2.1.vcf.bed|awk '{if($1==\"chr1\")print \"1A\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}' |awk '{if($1==\"chr2\")print \"1A\"\"\\t\"$2+471304005\"\\t\"$3+471304005\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr3\")print \"1B\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8; else print $0}'|awk '{if($1==\"chr4\")print \"1B\"\"\\t\"$2+438720154\"\\t\"$3+438720154\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr5\")print \"1D\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr6\")print \"1D\"\"\\t\"$2+452179604\"\\t\"$3+452179604\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr7\")print \"2A\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}' |awk '{if($1==\"chr8\")print \"2A\"\"\\t\"$2+462376173\"\\t\"$3+462376173\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr9\")print \"2B\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8; else print $0}'|awk '{if($1==\"chr10\")print \"2B\"\"\\t\"$2+453218924\"\\t\"$3+453218924\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr11\")print \"2D\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr12\")print \"2D\"\"\\t\"$2+462216879\"\\t\"$3+462216879\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr13\")print \"3A\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}' |awk '{if($1==\"chr14\")print \"3A\"\"\\t\"$2+454103970\"\\t\"$3+454103970\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr15\")print \"3B\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8; else print $0}'|awk '{if($1==\"chr16\")print \"3B\"\"\\t\"$2+448155269\"\\t\"$3+448155269\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr17\")print \"3D\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr18\")print \"3D\"\"\\t\"$2+476235359\"\\t\"$3+476235359\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr19\")print \"4A\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}' |awk '{if($1==\"chr20\")print \"4A\"\"\\t\"$2+452555092\"\\t\"$3+452555092\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr21\")print \"4B\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8; else print $0}'|awk '{if($1==\"chr22\")print \"4B\"\"\\t\"$2+451014251\"\\t\"$3+451014251\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr23\")print \"4D\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr24\")print \"4D\"\"\\t\"$2+451004620\"\\t\"$3+451004620\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr25\")print \"5A\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}' |awk '{if($1==\"chr26\")print \"5A\"\"\\t\"$2+453230519\"\\t\"$3+453230519\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr27\")print \"5B\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8; else print $0}'|awk '{if($1==\"chr28\")print \"5B\"\"\\t\"$2+451372872\"\\t\"$3+451372872\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr29\")print \"5D\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr30\")print \"5D\"\"\\t\"$2+451901030\"\\t\"$3+451901030\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr31\")print \"6A\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}' |awk '{if($1==\"chr32\")print \"6A\"\"\\t\"$2+452440856\"\\t\"$3+452440856\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr33\")print \"6B\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8; else print $0}'|awk '{if($1==\"chr34\")print \"6B\"\"\\t\"$2+452077197\"\\t\"$3+452077197\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr35\")print \"6D\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr36\")print \"6D\"\"\\t\"$2+450509124\"\\t\"$3+450509124\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr37\")print \"7A\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}' |awk '{if($1==\"chr38\")print \"7A\"\"\\t\"$2+450046986\"\\t\"$3+450046986\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr39\")print \"7B\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8; else print $0}'|awk '{if($1==\"chr40\")print \"7B\"\"\\t\"$2+453822637\"\\t\"$3+453822637\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr41\")print \"7D\"\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}'|awk '{if($1==\"chr42\")print \"7D\"\"\\t\"$2+453812268\"\\t\"$3+453812268\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8;else print $0}' >vmap2.1.vcf.bed");
                
//                    bw.newLine();
                      
//                bw.flush();
//                bw.close();
            }

//        } catch (Exception e) {
//            System.out.println("Error in calling sample GERP!");
//            e.printStackTrace();

        }

        //   String querySpecies[]={"hospo","aetau","brdis","ertef","hovul","leper","orbar","orbra","orgla","orglu","orind","orjap","orruf","orniv","orpun","ormer","orlon","ortho","pahal","phedu","seita","sobic","trura","zemay","pavir"};    
//         String querySpecies[]={"traesA","traesD","traesB"};//eccru
//        // String tree="(((((((((traesA aetau)trura))(hospo hovul))brdis)phedu)(leper(orbra(orpun(ormer(orlon(orglu((orbar orgla)((orniv orind)(orruf orjap))))))))))(((ertef(ercur ortho))((sobic zemay)(seita(eccru(pahal pavir)))))))ancom)";
//        String []subGenomes={"A","B","D"};
//        for(String subGenome:subGenomes){
//     //   for(String querySpecie:querySpecies){
//   
////    System.out.println("maf-swap /data2/lipeng/msa/"+querySpecie+"/traes"+subGenome+"."+querySpecie+".1split.maf >/data2/lipeng/msa/"+querySpecie+"/"+subGenome+".swap.maf");
////    System.out.println("java -Xms50g -Xmx50g -jar mafRevise.jar -m /data2/lipeng/msa/"+querySpecie+"/"+subGenome+".swap.maf -f /data1/home/lipeng/database/genome/"+querySpecie+"/genome.genome -o /data2/lipeng/msa/"+querySpecie+"/"+subGenome+".temp.maf -s "+querySpecie);
////    System.out.println("last-split -fMAF -m1 /data2/lipeng/msa/"+querySpecie+"/"+subGenome+".temp.maf |maf-swap >/data2/lipeng/msa/"+querySpecie+"/traes"+subGenome+"."+querySpecie+".2split.maf");
////   System.out.println("last-postmask "+"/data2/lipeng/msa/"+querySpecie+"/traes"+subGenome+"."+querySpecie+".2split.maf >/data2/lipeng/msa/"+querySpecie+"/traes"+subGenome+"."+querySpecie+".post.maf");
////    System.out.println("rm "+"/data2/lipeng/msa/"+querySpecie+"/"+subGenome+".temp.maf");
////    System.out.println("rm "+"/data2/lipeng/msa/"+querySpecie+"/"+subGenome+".swap.maf");
////    System.out.println("rm "+"/data2/lipeng/msa/"+querySpecie+"/traes"+subGenome+"."+querySpecie+".2split.maf");
////////    
////        System.out.println("cp ../"+querySpecie+"/traes"+subGenome+"."+querySpecie+".post.maf x.maf");
////         System.out.println("sed '/^#/d' x.maf >xx.maf");
////         System.out.println("cat head xx.maf >"+"traes"+subGenome+"."+querySpecie+".sing.maf");
////         System.out.println("rm xx.maf");
////          System.out.println("rm x.maf");
//  //      System.out.println("nohup roast E=traes"+subGenome+" \""+tree+"\""+" *.sing.maf 28way.maf");
////
////
////         }
// int []chrs={1,2,3,4,5,6,7};
////     System.out.println("maf-sort 28way."+subGenome+".maf >28way."+subGenome+".sort.maf");
//
//    for(int i:chrs){
//        System.out.println("mafSpeciesSubset "+i+subGenome+".sort.maf species.lst "+i+subGenome+".maf");
//        System.out.println("rm "+i+subGenome+".sort.maf");
//System.out.println("mafFilter -minRow=2 "+i+subGenome+".maf >"+i+subGenome+".flt.maf");
//System.out.println("java -Xms5g -Xmx5g -jar PlantGenetics.jar -m /data2/lipeng/asGenome/"+subGenome+"/"+i+subGenome+".flt.maf -o /data2/lipeng/asGenome/"+subGenome+"/"+i+subGenome+".maf");
//System.out.println("sed -n '"+i+"p' /data1/home/lipeng/database/genome/traes"+subGenome+"/genome.bed >/data1/home/lipeng/database/genome/traes"+subGenome+"/"+i+subGenome+".bed");
//System.out.println("bedtools getfasta -bed /data1/home/lipeng/database/genome/traes"+subGenome+"/"+i+subGenome+".bed -fi /data1/home/lipeng/database/genome/traes"+subGenome+"/genome.fa -fo /data1/home/lipeng/database/genome/traes"+subGenome+"/"+i+subGenome+".fa");
//System.out.println("prequel -k -r /data1/home/lipeng/database/genome/traes"+subGenome+"/"+i+subGenome+".fa -i MAF "+i+subGenome+".maf "+subGenome+".mod "+i+subGenome);
//System.out.println("prequel -k -n  -r /data1/home/lipeng/database/genome/traes"+subGenome+"/"+i+subGenome+".fa -i MAF "+i+subGenome+".maf "+subGenome+".mod "+i+subGenome);
//System.out.println("mafRanges /data2/lipeng/msa/gerp/"+i+subGenome+".sort.maf traes"+subGenome+" /data2/lipeng/msa/gerp/"+i+subGenome+".bed");
//System.out.println("java -jar gerpReformat.jar -f /data2/lipeng/msa/gerp/"+i+subGenome+".bed -m /data2/lipeng/msa/gerp/"+i+subGenome+".sort.maf.gerp -o /data2/lipeng/msa/gerp/"+i+subGenome+".gerp.wig");
//System.out.println("mv "+i+subGenome+".gerp.wigs "+i+subGenome+".gerp.wig");
//     System.out.println("sed -n '"+i+"p' "+subGenome+".bed >"+i+subGenome+".bed");
//     System.out.println("mafsInRegion "+i+subGenome+".bed "+i+subGenome+".maf 28way."+subGenome+".sort.maf");
//      System.out.println("maf-sort "+i+subGenome+".maf >" +i+subGenome+".sort.maf");
//      System.out.println("gerpcol -e traes"+j+" -z -t "+j+".gerp.tree -f " +i+j+".sort.maf"+" -x \".gerp\"");
//System.out.println("mafFilter -minRow=3 "+i+j+".sort.maf >"+i+j+".sort.flt.maf");
//System.out.println("mafRanges "+i+j+".sort.flt.maf traes"+j+" "+i+j+".flt.bed");  
//System.out.println("bedtools merge -i "+i+j+".flt.bed >"+i+j+".flt.beds");  
//System.out.println("rm "+i+j+".flt.bed");  
//System.out.println("mv "+i+j+".flt.beds "+i+j+".flt.bed");  
//System.out.println("java -Xms400g -Xmx400g -jar gerpReformat.jar -f /data2/lipeng/msa/gerp/"+i+j+".flt.bed -m /data2/lipeng/msa/gerp/"+i+j+".sort.maf.gerp -o /data2/lipeng/msa/gerp/"+i+j+".gerp.wig -g "+j);
//System.out.println("java -Xms200g -Xmx200g -jar PlantGenetics.jar -g /data2/lipeng/msa/gerp/"+i+j+".gerp.wig");

//   }  // 
//}
//(((orind,orniv),(orruf,orjap)),(orbar,orgla))))))leper)))
//(((((phedu,(((hovul,hospo),(aetau,(traesA,trura))),brdis)),(orbra,((orpun,(ormer,(orlon,(orglu,(((orind,orniv),(orruf,orjap)),(orbar,orgla)))))),leper)))),(((sobic,zemay),(eccru,(seita,(pahal,pavir)))),((ercur,ertef),ortho))),ancom)
//((((phedu,(((hovul,hospo),(aetau,(traesA,trura))),brdis)),(orbra,((orpun,(ormer,(orlon,(orglu,(((orind,orniv),(orruf,orjap),(orbar,orgla)))))),leper)))),(((sobic,zemay):,(eccru:,(seita,(pahal,pavir)))),((ercur,ertef),ortho))),ancom);
//(((((((((traesA aetau)trura))(hospo hovul))brdis)phedu)(leper(orbra(orpun(ormer(orlon(orglu((orbar orgla)((orniv orind)(orruf orjap))))))))))(((ertef(ercur ortho))((sobic zemay)(seita(eccru(pahal pavir)))))))ancom)
