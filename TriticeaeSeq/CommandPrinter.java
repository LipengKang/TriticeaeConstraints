/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

/**
 *
 * @author kanglipeng
 */
public class CommandPrinter {                                                                                                                                           
    
     public static void main(String[] args) {
         String querySpecies[]={"hospo","aetau","brdis","ertef","hovul","leper","orbar","orbra","orgla","orglu","orind","orjap","orruf","orniv","orpun","ormer","orlon","ortho","pahal","phedu","seita","sobic","trura","zemay","pavir"};    
        // String querySpecies[]={"eccru"};//eccru
         String tree="(((((((((traesA aetau)trura))(hospo hovul))brdis)phedu)(leper(orbra(orpun(ormer(orlon(orglu((orbar orgla)((orniv orind)(orruf orjap))))))))))(((ertef(ercur ortho))((sobic zemay)(seita(eccru(pahal pavir)))))))ancom)";
//        for(String querySpecie:querySpecies){
    String []subGenome={"A","B","D"};
////    System.out.println("maf-swap /data2/lipeng/msa/"+querySpecie+"/traes"+subGenome+"."+querySpecie+".1split.maf >/data2/lipeng/msa/"+querySpecie+"/"+subGenome+".swap.maf");
////    System.out.println("java -Xms50g -Xmx50g -jar mafRevise.jar -m /data2/lipeng/msa/"+querySpecie+"/"+subGenome+".swap.maf -f /data1/home/lipeng/database/genome/"+querySpecie+"/genome.genome -o /data2/lipeng/msa/"+querySpecie+"/"+subGenome+".temp.maf -s "+querySpecie);
////    System.out.println("last-split -fMAF -m1 /data2/lipeng/msa/"+querySpecie+"/"+subGenome+".temp.maf |maf-swap >/data2/lipeng/msa/"+querySpecie+"/traes"+subGenome+"."+querySpecie+".2split.maf");
////    System.out.println("last-postmask "+"/data2/lipeng/msa/"+querySpecie+"/traes"+subGenome+"."+querySpecie+".2split.maf >/data2/lipeng/msa/"+querySpecie+"/traes"+subGenome+"."+querySpecie+".post.maf");
////    System.out.println("rm "+"/data2/lipeng/msa/"+querySpecie+"/"+subGenome+".temp.maf");
////    System.out.println("rm "+"/data2/lipeng/msa/"+querySpecie+"/"+subGenome+".swap.maf");
////    System.out.println("rm "+"/data2/lipeng/msa/"+querySpecie+"/traes"+subGenome+"."+querySpecie+".2split.maf");
////    
////         System.out.println("cp ../"+querySpecie+"/traes"+subGenome+"."+querySpecie+".post.maf x.maf");
////         System.out.println("sed '/^#/d' x.maf >xx.maf");
////          System.out.println("cat head xx.maf >"+"traes"+subGenome+"."+querySpecie+".sing.maf");
////          System.out.println("nohup roast E=traes"+subGenome+" \""+tree+"\""+" *.sing.maf 28way.maf");
//
//
//         }
     int []chrs={1,2,3,4,5,6,7};
//     System.out.println("maf-sort 28way."+subGenome+".maf >28way."+subGenome+".sort.maf");
for(String j: subGenome){
     for(int i:chrs){
//     System.out.println("sed -n '"+i+"p' "+subGenome+".bed >"+i+subGenome+".bed");
//     System.out.println("mafsInRegion "+i+subGenome+".bed "+i+subGenome+".maf 28way."+subGenome+".sort.maf");
//      System.out.println("maf-sort "+i+subGenome+".maf >" +i+subGenome+".sort.maf");
//      System.out.println("gerpcol -e traes"+j+" -z -t "+j+".gerp.tree -f " +i+j+".sort.maf"+" -x \".gerp\"");
//System.out.println("mafFilter -minRow=3 "+i+j+".sort.maf >"+i+j+".sort.flt.maf");
//System.out.println("mafRanges "+i+j+".sort.flt.maf traes"+j+" "+i+j+".flt.bed");  
//System.out.println("bedtools merge -i "+i+j+".flt.bed >"+i+j+".flt.beds");  
//System.out.println("rm "+i+j+".flt.bed");  
//System.out.println("mv "+i+j+".flt.beds "+i+j+".flt.bed");  
System.out.println("java -Xms400g -Xmx400g -jar gerpReformat.jar -f /data2/lipeng/msa/gerp/"+i+j+".flt.bed -m /data2/lipeng/msa/gerp/"+i+j+".sort.maf.gerp -o /data2/lipeng/msa/gerp/"+i+j+".gerp.wig -g "+j);

     }}
     }  
}
//(((orind,orniv),(orruf,orjap)),(orbar,orgla))))))leper)))
//(((((phedu,(((hovul,hospo),(aetau,(traesA,trura))),brdis)),(orbra,((orpun,(ormer,(orlon,(orglu,(((orind,orniv),(orruf,orjap)),(orbar,orgla)))))),leper)))),(((sobic,zemay),(eccru,(seita,(pahal,pavir)))),((ercur,ertef),ortho))),ancom)
//((((phedu,(((hovul,hospo),(aetau,(traesA,trura))),brdis)),(orbra,((orpun,(ormer,(orlon,(orglu,(((orind,orniv),(orruf,orjap),(orbar,orgla)))))),leper)))),(((sobic,zemay):,(eccru:,(seita,(pahal,pavir)))),((ercur,ertef),ortho))),ancom);
//(((((((((traesA aetau)trura))(hospo hovul))brdis)phedu)(leper(orbra(orpun(ormer(orlon(orglu((orbar orgla)((orniv orind)(orruf orjap))))))))))(((ertef(ercur ortho))((sobic zemay)(seita(eccru(pahal pavir)))))))ancom)