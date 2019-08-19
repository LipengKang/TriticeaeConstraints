        /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import java.io.BufferedReader;
import java.util.ArrayList;

public class ExtractTools{
    public ExtractTools(){}
    void ParseCodonSnp(String vcf){

        System.out.println("vcf profiling...");
        String temp = null;
        String[] tem = null;
        BufferedReader br;
        if (!vcf.endsWith(".gz")) {
            br = utils.IOUtils.getTextReader(vcf);
        } else {
            br = utils.IOUtils.getTextGzipReader(vcf);
        }
        
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) 
                    continue;
                    
                    
                    
                    
                }}catch (Exception e) {
            System.out.println("Errors in vcf profiling! please try again! Come on !");
                e.printStackTrace();
        } 
        
        
    }
    
    
    
}
   


    
    


