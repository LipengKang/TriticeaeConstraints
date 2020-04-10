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
public class Ccsprofiler {
String specie="davil";
String fieldNum[]={"3","5","6","7","9","13","14","28","29","31"};
    public Ccsprofiler() {
        this.hiCanu();
    }

    // ccs preparation and assembly of Triticeae by hiCanu
    public void hiCanu() {
        
        BufferedWriter bw = IOUtils.getTextWriter("/Users/kanglipeng/Desktop/1A/hiCanu.sh");
        try {
//            bw.write("bamtools filter -tag \"np:>2\" -in "+specie+".ccs.bam | bamtools filter -tag \"rq:>0.99\"| bamtools convert -format fasta -out "+specie+".ccs.fasta");bw.newLine();
//            bw.write("awk -F \"/\" '/>/ {$0=\">specieNum_cell_\"$2} 1' "+specie+".ccs.fa >"+specie+".fa");
  for (String i :fieldNum) {          
bw.write("bgzip "+i+"_1.fa");bw.newLine();
bw.write("bgzip "+i+"_2.fa");bw.newLine();}
  
            bw.flush();
            bw.close();

        } catch (Exception e) {
            System.out.println("Error in generating hiCanu.sh,please check code from Ccsprofiler!");
            e.printStackTrace();

        }
    }

    public static void main(String[] args) {
        new Ccsprofiler();
    }
}
