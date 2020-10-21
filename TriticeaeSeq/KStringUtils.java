/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import java.io.File;
import java.util.ArrayList;
import utils.*;
import java.util.List;
import static utils.PStringUtils.fastSplit;

/**
 *
 * @author kanglipeng
 */
public class KStringUtils extends PStringUtils {

    /**
     * Return a list of split String using Guava splitter, blank space
     *
     * @param line
     * @return
     */
     public static List<String>  fastSplitSpa(String line) {
         List<String> ls = fastSplit(line, " ");
        return ls;
    }

    /**
     * Return a list of split String using Guava splitter, delimiter tab
     *
     * @param line
     * @return
     */
    public static String[] fastSplitdel(String line) {
        List<String> ls = fastSplit(line, "\t");
        String tem[] = null;
        tem = ls.toArray(new String[ls.size()]);
        return tem;
    }

    /**
     * Return a list of split String using Guava splitter, comma
     *
     * @param line
     * @return
     */
     public static String[] fastSplitComma(String line) {
        List<String> ls = fastSplit(line, ",");
         String tem[] = null;
        tem = ls.toArray(new String[ls.size()]);
        return tem;
    }

    /**
     * Return a list of split String using Guava splitter, strigula
     *
     * @param line
     * @return
     */

    public static String[] fastSplitStrigula(String line) {
        List<String> ls = fastSplit(line, "-");
        String tem[] = null;
        tem = ls.toArray(new String[ls.size()]);
        return tem;
    }

    /**
     * Return a list of split String using Guava splitter, dot
     *
     * @param line
     * @return
     */

    public static String[] fastSplitDot(String line) {
        List<String> ls = fastSplit(line, ".");
        String tem[] = null;
        tem = ls.toArray(new String[ls.size()]);
        return tem;
    }

    /**
     * Return a list of split String using Guava splitter, semicolon
     *
     * @param line
     * @return
     */

    public static List<String> fastSplitSemicolon(String line) {
        List<String> ls = fastSplit(line, ";");
        return ls;
    }

    /**
     * Return a list of split String using Guava splitter, colon
     *
     * @param line
     * @return
     */

    public static String[] fastSplitColon(String line) {
        List<String> ls = fastSplit(line, ":");
        String tem[] = null;
        tem = ls.toArray(new String[ls.size()]);
        return tem;
    }

    /**
     * Return a list of split String form maf esemble format
     *
     * @param line
     * @return
     */
    public static String[] mafSplit(String line) {
        List<String> ls = KStringUtils.fastSplitSpa(line);
        String tem[] = null;
        List<String> lsN = new ArrayList<>();
        for (int i = 0; i < ls.size(); i++) {
            if (ls.get(i) != null && !ls.get(i).equals("")) {
                lsN.add(ls.get(i));
            }
        }
        tem = lsN.toArray(new String[lsN.size()]);
        return tem;

    }

    /**
     * Return a un-empty Array of split String
     *
     * @param line
     * @return
     */

    public static void quickSort(int[] array, int left, int right) {
        if (left > right) {
            return;
        }
        // base中存放基准数
        int base = array[left];
        int i = left, j = right;
        while (i != j) {
            // 顺序很重要，先从右边开始往左找，直到找到比base值小的数
            while (array[j] >= base && i < j) {
                j--;
            }

            // 再从左往右边找，直到找到比base值大的数
            while (array[i] <= base && i < j) {
                i++;
            }

            // 上面的循环结束表示找到了位置或者(i>=j)了，交换两个数在数组中的位置
            if (i < j) {
                int tmp = array[i];
                array[i] = array[j];
                array[j] = tmp;
            }
        }

        // 将基准数放到中间的位置（基准数归位）
        array[left] = array[i];
        array[i] = base;
        // 递归，继续向基准的左右两边执行和上面同样的操作
        // i的索引处为上面已确定好的基准值的位置，无需再处理
        quickSort(array, left, i - 1);
        quickSort(array, i + 1, right);
    }

    /**
     * Return a number representing base ACGTN-> byte 01234
     *
     * @param base
     * @return
     */
    public static byte baseToByte(String base) {
        byte baseNum = 0;
        if (base.equals("A")||base.equals("a") ) {
            baseNum = 0;
        }
        if (base.equals("C")||base.equals("c")) {
            baseNum = 1;
        }
        if (base.equals("G")||base.equals("g")) {
            baseNum = 2;
        }
        if (base.equals("T")||base.equals("t")) {
            baseNum = 3;
        }
        if (base.equals("N")||base.equals("n")) {
            baseNum = 4;
        }
        return baseNum;
    }

    /**
     * Return files absolute path in path
     *
     * @param path
     * @return
     */
    public static String[] getFile(String path) {
        File file = new File(path);
        File[] array = file.listFiles();
        String[] filesPath = null;
        for (int i = 0; i < array.length; i++) {
            if (array[i].isFile()) {
                filesPath[i] = array[i].getPath();
            }
        }
        return filesPath;
    }
    
//    
//        public static int fastContain(List indexSet1,List set2,int position) {
//            
//    currentPos = Integer.parseInt(l.get(1));
//posIndex = startLists[chrIndex].binarySearch(currentPos);
//if (posIndex < 0) {
//posIndex = -posIndex-2;
//}
//if (posIndex < 0) continue;
//if (currentPos >= endLists[chrIndex].get(posIndex)) continue;}

}
