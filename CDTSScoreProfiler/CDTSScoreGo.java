/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CDTSScoreProfiler;

import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashLongIntMap;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import utils.CLIInterface;

/**
 *
 * @author lipeng
 * 
 */

public class CDTSScoreGo implements CLIInterface {
    Options options = new Options();
    HelpFormatter optionFormat = new HelpFormatter();
    String introduction = null;
    String mode = null;
    String kmerLengthS = null;
    String windowSizeS=null;
    String stropSizeS=null;
    String referenceGenomeFileS = null;

    String libFileS = null;
    String outputDirS = null;
    String vcfFile =null;
    int kmerLength = 0;
    int windowSize=0;
    int stropSize=0;
      HashIntIntMap[] intKmerMaps = null;
    HashLongIntMap[] longKmerMaps = null;
    
    public CDTSScoreGo (String[] args) {
        introduction = this.createIntroduction();
        this.createOptions();
        this.retrieveParameters (args);
        this.runProfiler();
    }

    void runProfiler () {
        if (mode.equals("b")) {
            ReferenceKmerLib lib = new ReferenceKmerLib(kmerLength, referenceGenomeFileS);
            lib.writeBinaryFile(libFileS);
        }
        else if (mode.equals("p")) {
            ReferenceKmerLib lib = new ReferenceKmerLib(kmerLength, referenceGenomeFileS);
            lib.writeBinaryFile(libFileS);
            
        if (kmerLength <= 16) {

            this.intKmerMaps = lib.intKmerMaps;
        } else {

            this.longKmerMaps = lib.longKmerMaps;
        }
            new GenomeProfiler(libFileS, referenceGenomeFileS,outputDirS, vcfFile,windowSize,stropSize,intKmerMaps,longKmerMaps);
            
        }
    }
    
    @Override
    public void retrieveParameters (String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            mode = line.getOptionValue("m");
            kmerLengthS = line.getOptionValue("k");
             windowSizeS = line.getOptionValue("w");
               stropSizeS = line.getOptionValue("s");
            referenceGenomeFileS = line.getOptionValue("r");
            libFileS = line.getOptionValue("l");
            vcfFile=line.getOptionValue("v");
            outputDirS = line.getOptionValue("o");
            
        }
        catch(Exception e) {
            e.printStackTrace();
        }
        if (mode == null) {
            this.printIntroductionAndUsage();
            return;
        }
        else if (mode.equals("b")) {
            if (kmerLengthS==null) {
                this.printIntroductionAndUsage();
                return;
            }
            else  {
                kmerLength = Integer.valueOf(kmerLengthS);
                
                return;
            }
        }
        else if (mode.equals("p")) {
            if (kmerLengthS == null&&stropSizeS==null&&windowSizeS==null) {
                this.printIntroductionAndUsage();
                return;
            }
            else {
                kmerLength = Integer.valueOf(kmerLengthS);
                stropSize=Integer.valueOf(stropSizeS);
                windowSize=Integer.valueOf(windowSizeS);
                return;
            }
        }
        else {
            this.printIntroductionAndUsage();
            return;
        }
    }
    
    @Override
    public void printIntroductionAndUsage () {
        System.out.println("Incorrect parameter input. Program quits.");
        System.out.println(introduction);
        optionFormat.printHelp("CDTSScoreProfiler.jar", options );
    }
    
    @Override
    public void createOptions () {
        options = new Options();
        options.addOption("m", true, "Analysis mode. Two modes are available, building reference kmer library (b option) and profiling CDTS score (p option). e.g. -m b");
        options.addOption("k", true, "Kmer length. Only 32 and 16 are supported. e.g. -k 32");
        options.addOption("w", true, "Constraint elements window size when scan CDTS score. e.g. -w 550");
        options.addOption("s", true, "window slide every stropSize bp. e.g. -s 10");
        options.addOption("r", true, "Reference genome file. e.g -i maizeAGPV4.fa");
        options.addOption("l", true, "Kmer library file. e.g -l maize_32mer.lib");
          options.addOption("v", true, "vcf file. e.g -v maizeAGPV4.vcf");
        options.addOption("o", true, "Output directory. e.g -o maizeAGPV4_CDTS");
    }
    
    public String createIntroduction () {
        StringBuilder sb = new StringBuilder();
        sb.append("\nThe program CDTSScoreProfiler.jar is designed to construct a constraint map using a given genome. ");
        sb.append("It has 2 analysis modes. The first is to build kmer library from the reference genome. The second is to calculate the CDTS Score from vcf.\n\n");
        sb.append("Command line example:\n\n");
        sb.append("\t1. Build kmer library from reference genome. The reference genome should be in Fasta format. The header must be the chromosome number. e.g. chromosome 1 is >1\n");
        sb.append("\t\tjava -jar CDTSScoreProfiler.jar -m b -k 32 -r maizeAGPV4.fa -l maize_32mer.lib\n");
        sb.append("\t2. Calculate AFScore for each kmer from vcf\n");
        sb.append("\t\tjava -jar CDTSScoreProfiler.jar -m p -k 32 -w 10 -s 5 -r maizeAGPV4.fa -l maize_32mer.lib -v maizeAGPV4.vcf -o maizeAGPV4_CDTS\n");
        return sb.toString();
    }
     public static void main(String[] args) {
        new CDTSScoreGo (args);
    }
}
