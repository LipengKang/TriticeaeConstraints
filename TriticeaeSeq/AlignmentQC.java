/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author kanglipeng map or alignment QC
 */
public class AlignmentQC {

    Options options = new Options();
    String gffDir = null;
    String mafDir = null;
    String subGenome = null;
    String annotation = null;
     String function = null;
    String temp = null;
    String tem[] = null;

    public AlignmentQC(String[] args) {   
        this.createOptions();
        this.retrieveParameters(args);
      //  this.createIntroduction();
        this.parseDensity(gffDir);

           //  this.diffMAF(gffDir, mafDir, subGenome);
         
    }

    //-----------------------------test maf cover rate of gene or exon or other annotations  ---------------------------------
    public void parseDensity(String gffDir) {
        int rangeMax=5;int rangeMin=-5;int binNum=20;
        BufferedReader brg = IOUtils.getTextReader(gffDir);
        ArrayList<Float>data=new ArrayList<>();
          try {
       while ((temp = brg.readLine()) != null) {
           if(temp.startsWith("f")) continue;
           data.add(Float.parseFloat(temp));
       }
       Float []dataFrame = data.toArray(new Float[data.size()]);
        float []frame=new float[data.size()];
         for(int i = 0; i < dataFrame.length; i++){
    frame[i] = dataFrame[i].floatValue();
}
      double[][]DPMatrix= CsvProvider.statsFloatPD(frame, rangeMin,rangeMax,binNum);
       System.out.println(gffDir);
      for(int i=0;i<binNum;i++){
       System.out.println(DPMatrix[i][0]+"\t"+DPMatrix[i][1]);
      }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void diffMAF(String gffDir, String mafDir, String subGenome) {
        BufferedReader brg = IOUtils.getTextReader(gffDir);//gffDir points to mafDir now 
        BufferedReader brm = IOUtils.getTextReader(mafDir);
        int[][] hitsMap = new int[7][850000000];
        int[][] hitsMap1 = new int[7][850000000];
        int hitsLength = 0;
        int chrIndex = 0;
        try {
            while ((temp = brg.readLine()) != null) {
                if (!temp.startsWith("s tr")) {
                    continue;
                }
                tem = KStringUtils.mafSplit(temp);
               // if (tem[1].endsWith(subGenome)) {
                    for (int q = Integer.parseInt(tem[2]); q <= Integer.parseInt(tem[2]) + Integer.parseInt(tem[3]) - 1; q++) {
                        chrIndex = Integer.parseInt(String.valueOf(tem[1].charAt(10))) - 1;
                        hitsMap[chrIndex][q] = 1;
                //   }
                }
            }
            brg.close();
            long coverScore = 0;
            long CoverBps = 0;
            while ((temp = brm.readLine()) != null) {
                if (!temp.startsWith("s tr")) {
                    continue;
                }
                tem = KStringUtils.mafSplit(temp);
                //second maf may contain only defined subgenome
                for (int p = Integer.parseInt(tem[2]); p <= Integer.parseInt(tem[2]) + Integer.parseInt(tem[3]) - 1; p++) {
                    coverScore++;
                    chrIndex = Integer.parseInt(String.valueOf(tem[1].charAt(10))) - 1;
                      hitsMap1[chrIndex][p] = 3;
                    if (hitsMap[chrIndex][p] == 1) {
                        hitsMap[chrIndex][p] = 2;
                    }

                }
            }
            brm.close();
             for (int j = 0; j < 7; j++) {
                 int hitsLength1=0;
                for (int i = 0; i < hitsMap1[j].length; i++) {
                    if (hitsMap1[j][i] == 3) {  
                        hitsLength1++;
                    } 
                }
                System.out.println(hitsLength1);
            }

            for (int j = 0; j < 7; j++) {
                for (int i = 0; i < hitsMap[j].length; i++) {
                    if (hitsMap[j][i] == 2) {
                        CoverBps++;
                        hitsLength++;
                    } else if (hitsMap[j][i] == 1) {
                        hitsLength++;
                    }
                }
                System.out.println(CoverBps);
                
            }
            DecimalFormat df = new DecimalFormat("0.0000");
            System.out.println("maf1:"+mafDir);
            System.out.println("maf2:"+gffDir);
            System.out.println("Ref maf covers " + hitsLength + " bp " + subGenome + " genome.");
            System.out.println("last maf covers " + coverScore + " bp " + subGenome + " genome.");
            System.out.println("overlapped block covers " + CoverBps + " bp  region.(" + 100 * Float.parseFloat(df.format((float) CoverBps / hitsLength)) + "%)");

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void createOptions() {
        options = new Options();
        options.addOption("m", true, "maf path");
        options.addOption("g", true, "gff3 Path");
        options.addOption("a", true, "annotaion type of checking. eg:check exon or gene aligned rates");
        options.addOption("s", true, "sub-Genome. eg:A");
        options.addOption("f", true, "function. eg:diffMAF");
    }

    public void retrieveParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            mafDir = line.getOptionValue("m");
            gffDir = line.getOptionValue("g");
            annotation = line.getOptionValue("a");
            subGenome = line.getOptionValue("s");
            function =line.getOptionValue("f");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void createIntroduction() {
        StringBuilder sb = new StringBuilder();
        sb.append("The program AlignmentQC.jar is designed to calculate alignment cover rate.\n");
        sb.append("\nCommand line example:\n");
        sb.append("getCoverRate: java -Xms10g -Xmx20g -jar AlignmentQC.jar -f getCoverRate -m test.maf -g test.gff3 -a exon -s A\n");
         sb.append("diffMaf: java -Xms10g -Xmx20g -jar AlignmentQC.jar -f diffMaf -m last.maf -g Ref.maf -s A\n");
        sb.append("\nAttention! " + "\n" + "1.maf must be one to one alignment (2-split or single_cov2 or chainNet) \n2.gff3 must contains only the longest transcripts(not affect gene cover rate) ");
        sb.append("\nResult of " + mafDir + ":\n");
        System.out.print(sb.toString());
    }

    public static void main(String[] args) {
        new AlignmentQC(args);
    }
}
