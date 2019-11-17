/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CDTSScoreProfiler;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashIntFloatMap;
import com.koloboke.collect.map.hash.HashIntFloatMaps;
import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import com.koloboke.collect.map.hash.HashLongFloatMap;
import com.koloboke.collect.map.hash.HashLongFloatMaps;
import com.koloboke.collect.map.hash.HashLongIntMap;
import com.koloboke.collect.map.hash.HashLongIntMaps;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import format.dna.BaseEncoder;
import format.dna.FastaByte;
import gnu.trove.list.array.TFloatArrayList;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
//import utils.Benchmark;
import utils.PArrayUtils;
import utils.PStringUtils;
import utils.IOUtils;
import CDTSScoreProfiler.CDTSScoreGo;
import java.text.DecimalFormat;
import CDTSScoreProfiler.VcfProfiler;

/**
 *
 * @author feilu & lipeng
 */
public class GenomeProfiler {

    int kmerLength = -1;
    int barcodeLength = -1;
    HashLongFloatMap[] longMaps = null;
    HashIntFloatMap[] intMaps = null;
//    HashLongIntMap[] longCountMaps = null;
    HashIntIntMap[] intCountMaps = null;
    HashIntFloatMap[] AFScoreMaps = null;
    HashIntFloatMap[] expectedCEScoreMaps = null;
    HashIntFloatMap[] CDTS = null;
    HashIntIntMap[] intKmerMaps = null;
    HashLongIntMap[] longKmerMaps = null;
    int stropSize = -1;
    int windowSize = -1;
    float fCDTSScore = 0;
    double dCDTSScore = 0;

    public GenomeProfiler(String libFileS, String referenceGenomeFileS, String outputDirS, String vcfFile, int windowSize, int stropSize,   HashIntIntMap[] intKmerMaps ,
    HashLongIntMap[] longKmerMaps) {
        
        FastaByte f = new FastaByte(referenceGenomeFileS);
     //   HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIByteMap();
        this.initializeKmerCountMap(libFileS);
        this.countKmers(f, vcfFile);
        this.callExpectedCEs(f,intKmerMaps,longKmerMaps);
        this.callCDTS(f, outputDirS, windowSize, stropSize);

    }

    public HashIntFloatMap[] callExpectedCEs(FastaByte f,HashIntIntMap[] intKmerMaps ,
    HashLongIntMap[] longKmerMaps) {
//        File outputDir = new File(outputDirS);
//        outputDir.mkdir();
        HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIByteMap();

        HashMap<Integer, String> chrSeqMap = new HashMap();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            chrSeqMap.put(Integer.valueOf(f.getName(i)), f.getSeq(i));
        }
        Set<Map.Entry<Integer, String>> chrSeqset = chrSeqMap.entrySet();
        System.out.println("calculate expected CDTS Score by chromosomes...");
        long start = System.nanoTime();
        this.expectedCEScoreMaps = new HashIntFloatMap[f.getSeqNumber()];
//        for (int i = 0; i < f.getSeqNumber(); i++) {
//            //make sure all chr have enough disk
//           // expectedCEScoreMaps[i] = HashIntFloatMaps.newMutableMap((int) (2 * f.getTotalSeqLength() / f.getSeqNumber()));
//           //in wheat lulab gff, all sub-chrs are smaller than 
//           expectedCEScoreMaps[i] = HashIntFloatMaps.newMutableMap((int) (f.);
//        }

        chrSeqset.parallelStream().forEach(entry -> {

            int chr = entry.getKey()-1;
           System.out.println("chromosome "+String.valueOf(chr+1)+" eCDTS calculating....");
            expectedCEScoreMaps[chr] = HashIntFloatMaps.newMutableMap((int) (f.getSeqLength(chr)));
            String seq = entry.getValue();
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
//            StringBuilder sb = new StringBuilder();
//            sb.append("chr").append(PStringUtils.getNDigitNumber(3, chr)).append("_AFScore.txt.gz");
//            String outCDTSfileS = new File(outputDir, sb2.toString()).getAbsolutePath();
//            BufferedWriter bw2 = IOUtils.getTextGzipWriter(outCDTSfileS);
//            try {
//
//                bw2.write("AFScore");
//                bw2.newLine();
            for (int i = 0; i < kmerLength / 2 - 1; i++) {
                expectedCEScoreMaps[chr].put(i, 1);
//                    bw2.write("NA");
//                    bw2.newLine();
            }
            int mark = 0;
            boolean flag = false;

            int stropMark = kmerLength / 2 - 1;
            if (intMaps != null) {

                for (int j = 0; j < bArray.length - kmerLength + 1; j++) {
                    flag = false;
                    for (int k = mark; k < j + kmerLength; k++) {
                        if (bArray[k] > 3) {
                            j = k;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        mark = j + 1;
//                            bw2.write("NA");
//                            bw2.newLine();
//                            bw2.write("NA");
//                            bw2.newLine();
                        continue;
                    } else {
                        mark = j + kmerLength;
                    }
                    int barcodeIndex = ReferenceKmerLib.getBarcodeIndex(bArray, j, barcodeLength);
                    int kmerV = BaseEncoder.getIntSeqFromSubByteArray(bArray, j, j + kmerLength);
                    if (intMaps[barcodeIndex].containsKey(kmerV)) {

                        float aver = intMaps[barcodeIndex].get(kmerV) / intKmerMaps[barcodeIndex].get(kmerV);
//                            bw2.write(String.valueOf(aver));
//                            bw2.newLine();
                        //initial all pos with N in adhere kmer
                         expectedCEScoreMaps[chr].put(j + kmerLength / 2 - 1, aver);
                        if (stropMark + 1 < j + kmerLength / 2 - 1) {
                            for (int i = stropMark; i < j + kmerLength / 2 - 1; i++) {
                                expectedCEScoreMaps[chr].put(i, 1);
                            }

                        }
                       

                        stropMark = j + kmerLength / 2 - 1;
                    }
                }

            } else if (longMaps != null) {

                for (int j = 0; j < bArray.length - kmerLength + 1; j++) {
                    flag = false;

                    for (int k = mark; k < j + kmerLength; k++) {
                        //define slidewindow with N in variants pos
                        if (bArray[k] > 3) {
                            j = k;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        mark = j + 1;
//                            bw2.write("NA");
//                            bw2.newLine();
//                            bw2.write("NA");
//                            bw2.newLine();

                        continue;
                    } else {
                        //no use to scan every pos for "N",scan new pos is enough
                        mark = j + kmerLength;
                    }
                    //barcode to find bin and then check kmerv whether or not exist
                    int barcodeIndex = ReferenceKmerLib.getBarcodeIndex(bArray, j, barcodeLength);
                    long kmerV = BaseEncoder.getLongSeqFromSubByteArray(bArray, j, j + kmerLength);
                    if (longMaps[barcodeIndex].containsKey(kmerV)) {

                        float aver = longMaps[barcodeIndex].get(kmerV) / longKmerMaps[barcodeIndex].get(kmerV);
//                            bw2.write(String.valueOf(aver));
//                            bw2.newLine();
expectedCEScoreMaps[chr].put(j + kmerLength / 2 - 1, aver);
                        if (stropMark + 1 < j + kmerLength / 2 - 1) {
                            for (int i = stropMark; i < j + kmerLength / 2 - 1; i++) {
                                expectedCEScoreMaps[chr].put(i, 1);
                            }

                        }
                    
                        stropMark = j + kmerLength / 2 - 1;
                    }

                    //                sb2.append("Writing copy number score from ").append(referenceGenomeFileS).append(" is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
                }

            }
            for (int i = 0; i < kmerLength / 2; i++) {
//                    bw2.write("NA");
//                    bw2.newLine();
                expectedCEScoreMaps[chr].put(stropMark + i + 1, 1);
            }
//                bw2.flush();
//                bw2.close();
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
        });
        StringBuilder sb = new StringBuilder();
    //    sb.append("calculate expected CE scores from geneome").append(" is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
       sb.append("calculate expected CE scores from geneome").append(" is finished.  ") ;
        System.out.println(sb.toString());
        return expectedCEScoreMaps;
    }

    void callCDTS(FastaByte f, String outputDirS, int windowSize, int stropSize) {
        File outputDir = new File(outputDirS);
        outputDir.mkdir();
        //calculate observed contraint socres

        DecimalFormat df = new DecimalFormat("0.000");
        this.CDTS = new HashIntFloatMap[f.getSeqNumber()];
        for (int i = 0; i < f.getSeqNumber(); i++) {
            CDTS[i] = HashIntFloatMaps.newMutableMap((int) f.getTotalSeqLength() * 2 / windowSize / f.getSeqNumber());
        }
        HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIByteMap();
        HashMap<Integer, String> chrSeqMap = new HashMap();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            chrSeqMap.put(Integer.valueOf(f.getName(i)), f.getSeq(i));
        }
        Set<Map.Entry<Integer, String>> chrSeqset = chrSeqMap.entrySet();
        System.out.println("calculating CDTSScore by chromosomes...");

        chrSeqset.parallelStream().forEach(entry -> {
            int chr = entry.getKey()-1 ;
            String seq = entry.getValue();
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }

            ArrayList<Integer> AFScorePosList = new ArrayList(AFScoreMaps[chr].keySet());
            Collections.sort(AFScorePosList);

            for (int j = 0; j < bArray.length; j++) {

                if (Collections.binarySearch(AFScorePosList, j) >= 0) {
                    dCDTSScore = 1000 * Math.log((double) (1 + AFScoreMaps[chr].get(j) - expectedCEScoreMaps[chr].get(j)));
                    fCDTSScore = Float.parseFloat(df.format((float) dCDTSScore));
                    CDTS[chr].put(j, fCDTSScore);
                } else {
                  //  float x=expectedCEScoreMaps[chr].get(j);
                    dCDTSScore = 1000 * Math.log((double) (1 + 1 - expectedCEScoreMaps[chr].get(j)));
                    fCDTSScore = Float.parseFloat(df.format((float) dCDTSScore));
                    CDTS[chr].put(j, fCDTSScore);
                }
            }
            System.out.print("chr");
            System.out.print(String.valueOf(chr + 1));
            System.out.println(" CDTSScore calling completed, start writing CDTS score.");
            StringBuilder sb1 = new StringBuilder();
            sb1.append(String.valueOf(kmerLength)).append("mer_").append("chr").append(PStringUtils.getNDigitNumber(3, chr+1)).append("_CDTSScore.txt.gz");
            String outCDTSfileS = new File(outputDirS, sb1.toString()).getAbsolutePath();
            BufferedWriter bw1 = IOUtils.getTextGzipWriter(outCDTSfileS);

            StringBuilder sb2 = new StringBuilder();
            sb2.append(String.valueOf(windowSize)).append("ww_").append(String.valueOf(stropSize)).append("sp_").append(String.valueOf(kmerLength)).append("mer_chr").append(PStringUtils.getNDigitNumber(3, chr+1)).append("_CDTSScore.txt.gz");
            String outWWfileS = new File(outputDirS, sb2.toString()).getAbsolutePath();
            BufferedWriter bw2 = IOUtils.getTextGzipWriter(outWWfileS);
            int mark = 1;
            fCDTSScore = 0;
            int windowStartPos = 0;
            try {
                for (int j = 0; j < bArray.length; j++) {
                    bw1.write(String.valueOf(CDTS[chr].get(j)));
                    bw1.newLine();
                }
                bw1.flush();
                bw1.close();
                for (int j = 0; j < bArray.length; j++) {

                    if (mark % windowSize == 0) {
                        bw2.write(String.valueOf(fCDTSScore));
                        bw2.newLine();
                        fCDTSScore = 0;
                        j = windowStartPos + stropSize - 1;
                        windowStartPos = j;
                        mark = 1;
                    } else {
                        fCDTSScore = fCDTSScore + CDTS[chr].get(j);
                        mark++;
                    }
                }

                bw2.flush();
                bw2.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    void countKmers(FastaByte f, String vcfFile) {
//        referenceGenomeFileS = "/Users/feilu/Documents/analysisL/pipelineTest/cpScore/maize_chr12.fa";
// here scan the same genome to add CDTS scores
//        VcfProfiler v = new VcfProfiler(vcfFile, f);

        this.AFScoreMaps = CDTSScoreProfiler.VcfProfiler.calculateAF(vcfFile, f);
        HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIByteMap();
        System.out.println("Start counting kmers in genome.");
        long start = System.nanoTime();
        if (intMaps != null) {
            for (int i = 0; i < f.getSeqNumber(); i++) {
                ArrayList<Integer> AFScorePosList = new ArrayList(AFScoreMaps[i].keySet());
                Collections.sort(AFScorePosList);
                byte[] bArray = f.getSeq(i).getBytes();
                for (int j = 0; j < bArray.length; j++) {
                    bArray[j] = ascIIByteMap.get(bArray[j]);
                }
                int mark = 0;
                boolean flag = false;
                for (int j = 0; j < bArray.length - kmerLength + 1; j++) {
                    flag = false;
                    for (int k = mark; k < j + kmerLength; k++) {
                        if (bArray[k] > 3) {
                            j = k;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        mark = j + 1;
                        continue;
                    } else {
                        mark = j + kmerLength;
                    }
                    int barcodeIndex = ReferenceKmerLib.getBarcodeIndex(bArray, j, barcodeLength);
                    int kmerV = BaseEncoder.getIntSeqFromSubByteArray(bArray, j, j + kmerLength);
                    if (Collections.binarySearch(AFScorePosList, j + kmerLength / 2 - 1) >= 0) {
                        // BigDecimal AFScore = new BigDecimal(String.valueOf(AFScoreMaps[i].get(j + kmerLength / 2-1)));
                        // BigDecimal oldKmerAFScore = new BigDecimal(String.valueOf(intMaps[barcodeIndex].get(kmerV)));
                        //BigDecimal kmerAFScoreBD = AFScore.add(oldKmerAFScore);
                        intMaps[barcodeIndex].addValue(kmerV, AFScoreMaps[i].get(j + kmerLength / 2 - 1));
//                            kmerCounts = longCountMaps[barcodeIndex].get(kmerV) + 1;
                        //intMaps[barcodeIndex].put(kmerV, kmerAFScoreBD.floatValue());
//                            longCountMaps[barcodeIndex].addValue(kmerV, kmerCounts);
                    } else {
                        intMaps[barcodeIndex].addValue(kmerV, 1);
                    }

                    int pos = j + 1;
                    if (pos % 50000000 == 0) {
                        System.out.println("Chromosome: " + f.getName(i) + ". Length = " + String.valueOf(bArray.length) + "bp. Position: " + String.valueOf(pos));
                    }
                }
            }
        } else if (longMaps != null) {
            for (int i = 0; i < f.getSeqNumber(); i++) {
                ArrayList<Integer> AFScorePosList = new ArrayList(AFScoreMaps[i].keySet());
                Collections.sort(AFScorePosList);
                byte[] bArray = f.getSeq(i).getBytes();
                for (int j = 0; j < bArray.length; j++) {
                    bArray[j] = ascIIByteMap.get(bArray[j]);
                }
                int mark = 0;
                boolean flag = false;
                for (int j = 0; j < bArray.length - kmerLength + 1; j++) {
                    flag = false;

                    for (int k = mark; k < j + kmerLength; k++) {
                        if (bArray[k] > 3) {
                            j = k;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        mark = j + 1;
                        continue;
                    } else {
                        //no use to scan every pos for "N",scan new pos is enough
                        mark = j + kmerLength;
                    }
                    //barcode to find bin and then check kmerv whether or not exist
                    int barcodeIndex = ReferenceKmerLib.getBarcodeIndex(bArray, j, barcodeLength);
                    long kmerV = BaseEncoder.getLongSeqFromSubByteArray(bArray, j, j + kmerLength);
//                    int kmerCounts = 0;
                    if (longMaps[barcodeIndex].containsKey(kmerV)) {
                        //refactor here ! need extract vcf informations to calculate CDTS 
                        //32mer start point i ,snp should locate at i+15
                        // i-->chr,  k-->strop point, j-->0-based kmer start point
//                        kmerCounts = longCountMaps[barcodeIndex].get(kmerV) + 1;
//                        longCountMaps[44].addValue(kmerV, kmerCounts);
                        if (Collections.binarySearch(AFScorePosList, j + kmerLength / 2 - 1) >= 0) {
//                            BigDecimal AFScore = new BigDecimal(String.valueOf(AFScoreMaps[i].get(j + kmerLength / 2 - 1)));
//                            BigDecimal oldKmerAFScore = new BigDecimal(String.valueOf(longMaps[barcodeIndex].get(kmerV)));
//                            BigDecimal kmerAFScoreBD = AFScore.add(oldKmerAFScore);
//                            longMaps[barcodeIndex].put(kmerV, kmerAFScoreBD.floatValue());
                            longMaps[barcodeIndex].addValue(kmerV, AFScoreMaps[i].get(j + kmerLength / 2 - 1));
//                            longCountMaps[barcodeIndex].addValue(kmerV, kmerCounts);
                        } else {
                            longMaps[barcodeIndex].addValue(kmerV, 1);
//                            BigDecimal AFScore = new BigDecimal(String.valueOf(1));
//                            BigDecimal oldKmerAFScore = new BigDecimal(String.valueOf(longMaps[barcodeIndex].get(kmerV)));
//                            BigDecimal kmerAFScoreBD = AFScore.add(oldKmerAFScore);
//                            longMaps[barcodeIndex].put(kmerV, kmerAFScoreBD.floatValue());
                        }
                    }
                    //Reverse complement may not be suitable to CDTS
//                    long rKmerV = BaseEncoder.getLongReverseComplement(kmerV, kmerLength);
//                    if (longMaps[barcodeIndex].containsKey(rKmerV)) {
//                        longMaps[barcodeIndex].addValue(rKmerV, 1);
//                    }
                    int pos = j + 1;
                    if (pos % 50000000 == 0) {
                        System.out.println("Chromosome: " + f.getName(i) + ". Length = " + String.valueOf(bArray.length) + "bp. Position: " + String.valueOf(pos));
                    }
                }
            }
        }
        StringBuilder sb = new StringBuilder();
       // sb.append("Couting kmer is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
       sb.append("Couting kmer is finished.");
       System.out.println(sb.toString());
    }

    void initializeKmerCountMap(String libFileS) {
//        libFileS = "/Users/feilu/Documents/analysisL/pipelineTest/cpScore/kmerLib/kmerLib.bin"; 
        System.out.println("Reading kmer library from " + libFileS);
        long start = System.nanoTime();
        try {
            DataInputStream dis = IOUtils.getBinaryReader(libFileS);
//binary lib file  format       
//kmerLength,barcodeLength,kmerBin size,int or long kmers
            kmerLength = dis.readInt();
            barcodeLength = dis.readInt();
            int setSize = (int) (Math.pow(4, barcodeLength));
            if (kmerLength > 16) {
                longMaps = new HashLongFloatMap[setSize];
//                longCountMaps = new HashLongIntMap[setSize];
            } else if (kmerLength <= 16) {
                intMaps = new HashIntFloatMap[setSize];
//                intCountMaps = new HashIntIntMap[setSize];
            }
            int mapSize = -1;
            for (int i = 0; i < setSize; i++) {
                mapSize = dis.readInt();
                System.out.println("Map size of index " + String.valueOf(i) + ": " + mapSize);
                float[] values = new float[mapSize];
//                int[] countValues = new int[mapSize];
                if (kmerLength > 16) {
                    long[] keys = new long[mapSize];
                    for (int j = 0; j < mapSize; j++) {
                        keys[j] = dis.readLong();
                    }
                    longMaps[i] = HashLongFloatMaps.newMutableMap(keys, values);
//                    longCountMaps[i] = HashLongIntMaps.newImmutableMap(keys, countValues);
                } else if (kmerLength <= 16) {
                    int[] keys = new int[mapSize];
                    for (int j = 0; j < mapSize; j++) {
                        keys[j] = dis.readInt();
                    }
                    intMaps[i] = HashIntFloatMaps.newMutableMap(keys, values);
//                    intCountMaps[i]=HashIntIntMaps.newImmutableMap(keys, countValues);
                }
            }
            dis.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
       // System.out.println(Benchmark.getTimeSpanSeconds(start) + " seconds used to initialize the KmerCount map");
       System.out.println("complete initializing the KmerCount map");}
}
