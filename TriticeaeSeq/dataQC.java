/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import daxing.md5.MD5;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import utils.IOUtils;

/**
 *
 * @author kanglipeng
 */
//Triticeae PacBio Seq data QC
public class dataQC {

    int sum = 0;

    public dataQC() {
        this.checkMd5();

    }

//---------------------------------checkMd5----------------------------------------
//this part results a list of file need md5 and system out all checked boolean 
    public void checkMd5() {
        File md5txt = new File("/data1/home/lipeng/database/T.urartu_pacbio/sub_fastq/All_md5.txt");
        //create .txt to list all aboulute files path in a directory
       /* File file = new File(fileDir);
        File[] objectList = file.listFiles();
        ArrayList<String> fileList = new ArrayList();
        for (int i = 0; i < objectList.length; i++) {
            if (objectList[i].isFile()) {
                fileList.add(objectList[i].getAbsolutePath());
            }
        }

        //      new File(fileDir, "/md5list.txt").mkdirs();
        StringBuilder tmpFile = new StringBuilder();
        tmpFile.append(fileDir).append("/md5list.txt");

        try {
            BufferedWriter bw = null;
            bw = IOUtils.getTextWriter(tmpFile.toString());
            for (String i : fileList) {
                bw.write(i);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        File hashFileS = new File(tmpFile.toString());*/
        MD5.checkMD5ForDir(md5txt);
    }

    public static void main(String[] args) {
        new dataQC();
    }
}
