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
        try {
            //create blocks for all chr
            int chr=7;
            String genome="D";
            for (int j =1 ; j <= 64; j++) {
                String outDir = "/Users/kanglipeng/Desktop/1A/" + j + ".sh";                     
                BufferedWriter bw = IOUtils.getTextWriter(outDir);
                //block num
                int start = 2000*(j-1)+1;
                int end = j*2000;
                bw.write("for i in {" + start + ".." + end + "};do echo $i >$i.bed;done");
                bw.newLine();
                //question
                bw.write("for i in {" + start + ".." + end + "};do awk '{print \"chr"+chr+"\"\"\\t\"$1*5000-5000\"\\t\"$1*5000+10000}' $i.bed >$i.beds ;done");
                bw.newLine();
                
                bw.write("grep 'chr"+chr+"' ../"+genome+"/final/"+genome+".pns.bed  >"+chr+genome+".bed");bw.newLine();
                for (int i = start; i <= end-2; i++) {
                    //  System.out.println("bedtools intersect -a "+i+".bed -b 1A.neutral.bed > "+i+".block.bed");
//         System.out.println("awk '{if(($3-$2)>=100){print$0}}' "+i+".block.bed >"+i+".block.beds");
//         //find . -size 0 -delete 
                    //System.out.println("nohup phyloP --mode CONACC --wig-scores A.10k.mod ../gerp/"+i+"A.sort.maf >"+i+"A.phyloP.wig &");

                    bw.write("bedtools intersect -a " + i + ".beds -b ../"+genome+"/final/"+chr+genome+"/"+j+".bed  >" + i + ".bed");
                    bw.newLine();
                    // bw.write("mv "+i+".beds "+i+".bed");bw.newLine(); 
                   // bw.newLine();
                    bw.write("awk '{sum+=($3-$2)};END{print sum}' " + i + ".bed |awk '{if($1<100){cmd=\"rm " + i + ".bed\"; print cmd; system(cmd)}}'");
                    bw.newLine();
                   // bw.write("awk '{sum+=($3-$2)};END{print sum}' " + i + ".bed |awk '{if($1<100){cmd=\"rm " + i + ".beds\"; print cmd; system(cmd)}}'");
                   // bw.newLine();
                    //1.flt.maf

                    bw.write("mafsInRegion " + i + ".bed " + i + ".maf ../asGenome/"+genome+"/split/"+chr+genome+"/"+j+".maf");
                    bw.newLine();
                    //                  here need change
                    bw.write("phyloFit --tree \"(((traesD,aetau),trura),hovul)\" --subst-mod REV --out-root " + i + " " + i + ".maf");
                    bw.newLine();
                      bw.write("grep 'aetau' " + i + ".mod >"+i+".mods");
                                          bw.newLine();
                    //traesA and trura and 1A
                    // bw.write("all_dists -m "+i+".mod | awk '/traesA/{print$0}' | awk '/trura/{print$3}' >>1A.lemda");bw.newLine();
                    //                                                                                                                 here need change
                    bw.write("grep 'aetau' " + i + ".mods |sed 's/(/ /g'|sed 's/)/ /g'|sed 's/;/ /g'|sed 's/,/ /g'|xargs -n1|awk '/traesD/{print substr($1,8,length($1))}' >>" + j + ".lambda");
                    bw.newLine();
                    bw.write("rm " + i + ".beds");
                    bw.newLine();
                    bw.write("rm " + i + ".maf");
                    bw.newLine();
                      bw.write("rm " + i + ".bed");
                    bw.newLine();
                    bw.write("mv " + i + ".mods "+i+".mod");
                    bw.newLine();
                }
                      for (int i = end-1; i <= end; i++) {
                        bw.write("bedtools intersect -a " + i + ".beds -b "+chr+genome+".bed >" + i + ".bed");
                    bw.newLine();
                    bw.write("awk '{sum+=($3-$2)};END{print sum}' " + i + ".bed |awk '{if($1<100){cmd=\"rm " + i + ".bed\"; print cmd; system(cmd)}}'");
                    bw.newLine();
                    bw.write("mafsInRegion " + i + ".bed " + i + ".maf ../asGenome/"+genome+"/split/"+genome+".maf");
                    bw.newLine();
                    bw.write("phyloFit --tree \"(((traesD,aetau),trura),hovul)\" --subst-mod REV --out-root " + i + " " + i + ".maf");
                    bw.newLine();
                    
                         bw.write("grep 'aetau' " + i + ".mod >"+i+".mods");
                       bw.newLine();
                    
                    bw.write("grep 'aetau' " + i + ".mod |sed 's/(/ /g'|sed 's/)/ /g'|sed 's/;/ /g'|sed 's/,/ /g'|xargs -n1|awk '/traesD/{print substr($1,8,length($1))}' >>" + j + ".lambda");
                    bw.newLine();
                    bw.write("rm " + i + ".beds");
                    bw.newLine();
                    bw.write("rm " + i + ".maf");
                    bw.newLine();
                    bw.write("rm " + i + ".bed");
                    bw.newLine();
                    bw.write("mv " + i + ".mods "+i+".mod");
                    bw.newLine();
                      }
                bw.flush();
                bw.close();
            }

        } catch (Exception e) {
            System.out.println("Error in calling sample GERP!");
            e.printStackTrace();

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
    }
}
//   }  // 
//}
//(((orind,orniv),(orruf,orjap)),(orbar,orgla))))))leper)))
//(((((phedu,(((hovul,hospo),(aetau,(traesA,trura))),brdis)),(orbra,((orpun,(ormer,(orlon,(orglu,(((orind,orniv),(orruf,orjap)),(orbar,orgla)))))),leper)))),(((sobic,zemay),(eccru,(seita,(pahal,pavir)))),((ercur,ertef),ortho))),ancom)
//((((phedu,(((hovul,hospo),(aetau,(traesA,trura))),brdis)),(orbra,((orpun,(ormer,(orlon,(orglu,(((orind,orniv),(orruf,orjap),(orbar,orgla)))))),leper)))),(((sobic,zemay):,(eccru:,(seita,(pahal,pavir)))),((ercur,ertef),ortho))),ancom);
//(((((((((traesA aetau)trura))(hospo hovul))brdis)phedu)(leper(orbra(orpun(ormer(orlon(orglu((orbar orgla)((orniv orind)(orruf orjap))))))))))(((ertef(ercur ortho))((sobic zemay)(seita(eccru(pahal pavir)))))))ancom)
