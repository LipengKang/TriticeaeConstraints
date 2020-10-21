package CDTSScoreProfiler;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import com.koloboke.collect.map.hash.HashLongIntMap;
import com.koloboke.collect.map.hash.HashLongIntMaps;
import com.koloboke.collect.set.hash.HashIntSet;
import com.koloboke.collect.set.hash.HashIntSets;
import com.koloboke.collect.set.hash.HashLongSet;
import com.koloboke.collect.set.hash.HashLongSets;
import java.io.DataOutputStream;
import format.dna.BaseEncoder;
import format.dna.FastaByte;
//import utils.Benchmark;
import utils.IOUtils;

public class ReferenceKmerLib {

    int kmerLength = -1;
    HashIntSet[] intSets = null;
    HashLongIntMap[] longKmerMaps = null;
    HashIntIntMap[] intKmerMaps = null;
    int barcodeLength = 3;
    int binSize = (int) (Math.pow(4, barcodeLength));
 

    public ReferenceKmerLib(int kmerLength, String inputGenomeFileS) {
        this.kmerLength = kmerLength;
        this.createKmerSet(inputGenomeFileS);
    }

    void createKmerSet(String referenceGenomeFileS) {
        FastaByte f = new FastaByte(referenceGenomeFileS);
        HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIByteMap();
//getAscIIByteMap function convert ATCG to 0123,and others to 4，maybe aim to create new ATCG specific ASCII code 
        System.out.println("Building kmer library from reference...");
        System.out.println("KmerLength = " + String.valueOf(kmerLength) + " bp");
        long start = System.nanoTime();
        if (kmerLength <= 16) {
            this.intKmerMaps = this.getIntKmerMaps(f, ascIIByteMap);
        } else if (kmerLength > 16) {
            this.longKmerMaps = this.getLongKmerMaps(f, ascIIByteMap);
        }
 // System.out.println(Benchmark.getTimeSpanSeconds(start) + " seconds used to build Kmer library");
    }

    void writeBinaryFile(String libFileS) {
//kmer binary library file format       
//kmerLength, barcodeLength, kmerBin size,int(16) or long(32) kmers
        try {
            DataOutputStream dos = IOUtils.getBinaryWriter(libFileS);
            dos.writeInt(kmerLength);
            dos.writeInt(barcodeLength);
      
            if (kmerLength <= 16) {
                for (int i = 0; i < intKmerMaps.length; i++) {
                    int[] kmerD = intKmerMaps[i].keySet().toIntArray();
                    dos.writeInt(kmerD.length);
                    for (int j = 0; j < kmerD.length; j++) {
                        dos.writeInt(kmerD[j]);
                    }
                }
            } else if (kmerLength > 16) {
                for (int i = 0; i < longKmerMaps.length; i++) {
                    long[] kmerD = longKmerMaps[i].keySet().toLongArray();
                    dos.writeInt(kmerD.length);
                    for (int j = 0; j < kmerD.length; j++) {
                        dos.writeLong(kmerD[j]);
                    }
                }
            }
            dos.flush();
            dos.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Kmer library is written to " + libFileS);
    }

    public static int getBarcodeIndex(byte[] seqByte, int startIndex, int barcodeLength) {
        int index = 0;
        for (int i = 0; i < barcodeLength; i++) {
            int base = 0;
            if (i == 0) {
                base = seqByte[i + startIndex];
            } else {
                base = (int) (Math.pow(4, i)) * seqByte[i + startIndex];
            }
            index += base;
        }
        return index;
    }

    public  HashLongIntMap[] getLongKmerMaps(FastaByte f, HashByteByteMap ascIIByteMap) {

        int genomeSize = (int) f.getTotalSeqLength();
        HashLongIntMap[] longKmerMaps = new HashLongIntMap[binSize];
        for (int i = 0; i < binSize; i++) {
            longKmerMaps[i] = HashLongIntMaps.newMutableMap(genomeSize / binSize);
        }
        int barcodeIndex = 0;
        for (int k = 0; k < f.getSeqNumber(); k++) {
            String seq = f.getSeq(k);
            //sequence convert to byte further convert to former defined ASCII
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
            int mark = 0;
            boolean flag = false;
            for (int i = 0; i < bArray.length - kmerLength + 1; i++) {
                flag = false;
                //scan ambiguous letters---whenever "N" exists in 32 bps, move start point (i) to strop point (j)
                //question? such sequence like AAAAAAAAAAAAA......AAAANNNNNNNNNNAAAAAAAA.....      strop moves to N, i follows, marker move i's next step, this means all 32mers with "N" are discarded
                for (int j = mark; j < i + kmerLength; j++) {
                    if (bArray[j] > 3) {
                        i = j;
                        flag = true;
                        break;
                    }
                }

                if (flag) {
                    mark = i + 1;
                    continue;
                } else {
                    mark = i + kmerLength;
                }
                barcodeIndex = this.getBarcodeIndex(bArray, i, barcodeLength);
                long kmerL = BaseEncoder.getLongSeqFromSubByteArray(bArray, i, i + kmerLength);

                //64(4^3) bins as barcode to allocate kmers
                if (longKmerMaps[barcodeIndex].containsKey(kmerL)) {
                    longKmerMaps[barcodeIndex].addValue(kmerL, 1);
                } else {
                    longKmerMaps[barcodeIndex].put(kmerL, 1);
                }

                int pos = i + 1;
                //display kmer establish stats in every 50M windown
                if (pos % 50000000 == 0) {
                    long total = 0;
                    for (int u = 0; u < longKmerMaps.length; u++) {
                        total += (long) longKmerMaps[u].size();
                    }
                    System.out.println("Chromosome: " + f.getName(k) + ". Length = " + String.valueOf(bArray.length) + "bp. Position: " + String.valueOf(pos) + ". Kmer set size: " + String.valueOf(total));

                }
            }
        }
        return longKmerMaps;
    }

    public HashIntIntMap[] getIntKmerMaps(FastaByte f, HashByteByteMap ascIIByteMap) {
  
        int genomeSize = (int) f.getTotalSeqLength();
        HashIntIntMap[] intKmerMaps = new HashIntIntMap[binSize];
        //create (size=4^barcodeLength) array kmerSets 

        for (int i = 0; i < binSize; i++) {
            intKmerMaps[i] = HashIntIntMaps.newMutableMap(genomeSize / binSize);
        }
        //kmerSets[1]   kmerSets[2]  kmerSets[3]   ......  kemrSets[4^barcodeLength]
        //(.......)          
        //genomeSize/setSize
        int barcodeIndex = 0;
        for (int k = 0; k < f.getSeqNumber(); k++) {
            String seq = f.getSeq(k);
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
            int mark = 0;
            boolean flag = false;
            for (int i = 0; i < bArray.length - kmerLength + 1; i++) {
                flag = false;
                for (int j = mark; j < i + kmerLength; j++) {
                    if (bArray[j] > 3) {
                        i = j;
                        flag = true;
                        break;
                    }
                }
                if (flag) {
                    mark = i + 1;
                    continue;
                } else {
                    mark = i + kmerLength;
                }
                barcodeIndex = this.getBarcodeIndex(bArray, i, barcodeLength);
                int kmeri = BaseEncoder.getIntSeqFromSubByteArray(bArray, i, i + kmerLength);
                if (intKmerMaps[barcodeIndex].containsKey(kmeri)) {
                    intKmerMaps[barcodeIndex].addValue(kmeri, 1);
                } else {
                    intKmerMaps[barcodeIndex].put(kmeri, 1);
                }

                int pos = i + 1;
                if (pos % 50000000 == 0) {
                    long total = 0;
                    for (int u = 0; u < intKmerMaps.length; u++) {
                        total += (long) intKmerMaps[u].size();
                    }
                    System.out.println("Chromosome: " + f.getName(k) + ". Length = " + String.valueOf(bArray.length) + "bp. Position: " + String.valueOf(pos) + ". Kmer set size: " + String.valueOf(total));
                }
            }
        }
        return intKmerMaps;
    }
}
