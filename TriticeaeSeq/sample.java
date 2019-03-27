/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import static lipengKang.analysis.RandomArray.randomLongCommon;
import utils.PStringUtils;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author kanglipeng
 */
public class sample {

    String temp = null;
    String tem[] = null;

    public sample(String GERPDir, String outFile) {
        this.sampleGERP(GERPDir, outFile);
    }

    public void sampleGERP(String GERPDir, String outFile) {
        long wheatWholeGenome = 4934891648l;
int max=10000;
        long[] randomPos = randomLongCommon(1, wheatWholeGenome, max);
        Arrays.sort(randomPos);
        long pos = 1;
        int sampleIndex = 0;

        ArrayList<String> GERPsample = new ArrayList();
        BufferedReader br;
        BufferedWriter bw;
        if (GERPDir.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(GERPDir);
            bw = YaoIOUtils.getTextGzipWriter(outFile);
        } else {
            br = YaoIOUtils.getTextReader(GERPDir);
            bw = YaoIOUtils.getTextWriter(outFile);
        }
        try {
            while ((temp = br.readLine()) != null) {
                if (pos == randomPos[sampleIndex]) {
                    List<String> tList = PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    GERPsample.add(tem[1]);
                    sampleIndex++;
                    if (sampleIndex == max) {
                        break;
                    }
                }
                pos++;
            }
            for (String i : GERPsample) {
                bw.write(i);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in calling sample GERP!");
            e.printStackTrace();
            
        }

    }

    public static void main(String[] args) {
        String GERPDir = "/data1/home/lipeng/result/GERP/axt/D/gerp/wheatD.gerp++";
        String outFile = "/data1/home/lipeng/result/GERP/axt/D/gerp/wheatD.sample.gerp++";
        new sample(GERPDir, outFile);
    }

}
