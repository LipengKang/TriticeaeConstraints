/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TriticeaeSeq;

import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author kanglipeng
 */
public class CsvProvider {

    /**
     *  
     *  * statistic for probability density  @param frames:dataset 
     *  @param intervalCount:number of bins  @return 
     */
    public static double[][] statisticIntPD(int[] frames, int intervalCount) {
//probability density
        double Pd;
//bin parameters
        int upInterval, downInterval, middleValue;
// numbers in bin  
        int count = 0;
// 横纵坐标保存数组  
        double[][] frameArray = new double[intervalCount][2];

        Arrays.sort(frames);
        int framenum = frames.length;
        int minFrame = frames[0];
        int maxFrame = frames[framenum - 1];
// bin size  
        int interval = (maxFrame - minFrame) / intervalCount;
        System.out.println("Min=" + minFrame + " " + "Max=" + maxFrame);
        for (int k = 0; k < intervalCount; k++) {
            upInterval = minFrame + (k + 1) * interval;
            downInterval = minFrame + k * interval;
            middleValue = downInterval + interval / 2;
            for (int i = 0; i < framenum; i++) {
                if (frames[i] < upInterval && frames[i] >= downInterval) {
                    count++;
                }
            }
            Pd = (double) count / framenum / interval;
            frameArray[k][0] = middleValue;
            frameArray[k][1] = Pd;
            count = 0;
        }

        return frameArray;
    }

    /**
     *  
     *  * statistic for probability density 
     *  @param frames:dataset 
     *  @param range   @param binCount:number of bins  @return 
     */
    public static double[][] statsFloatPD(float[] frames, int min, int max, int binCount) {
//probability density
        double Pd;
//bin parameters
        float upInterval, downInterval;
        float middleValue;
// numbers in bin  
        int count = 0;
// 横纵坐标保存数组  
        double[][] frameArray = new double[binCount][2];

// bin size  
        float binSize = Float.valueOf(max - min) / binCount;
        float minData=0l;
            float maxData=0l;
        for (int k = 0; k < binCount; k++) {
            upInterval = min + (k + 1) * binSize;
            downInterval = min + k * binSize;
            middleValue = downInterval + binSize / 2;
            int frameNum = frames.length;
            
            for (int i = 0; i < frameNum; i++) {
                if(frames[i]>maxData){maxData=frames[i];}
                if(frames[i]<minData){minData=frames[i];}
                if (frames[i] < upInterval && frames[i] >= downInterval) {
                    count++;
                }
                if(k==0&&frames[i]<downInterval){count++;}
                if(k==binCount-1 && frames[i]>upInterval){count++;}
            }
            Pd = (double) count / frameNum ;
            BigDecimal b = new BigDecimal(Pd);
            double scalePd = b.setScale(2, BigDecimal.ROUND_HALF_UP).doubleValue();
            frameArray[k][0] = middleValue;
            frameArray[k][1] = scalePd;
            count = 0;
        }
   System.out.println("Min=" + minData + "\t" + "Max=" + maxData);
   System.out.println("Range=[" + min + "," + max+"]");
        return frameArray;
    }
    
    

/**
     *  
     *  statistic for standard deviation (sd)
     *  @param array :dataset
     * @return 
     */
	public static double statsDoubleSD(double[]array) {
		int sum = 0;      
		for(int i=0;i<array.length;i++){
		    sum += array[i];      //求出数组的总和
		}
		System.out.println(sum);  //939
		double average = sum/array.length;  //求出数组的平均数
		System.out.println(average);   //52.0
		int total=0;
		for(int i=0;i<array.length;i++){
		    total += (array[i]-average)*(array[i]-average);   //求出方差，如果要计算方差的话这一步就可以了
		}
		double standardDeviation = Math.sqrt(total/array.length);   //求出标准差
		return standardDeviation;
	}
        
 /**
     *  
     *  statistic for boxplot q1,q2,q3,q4,q5,q6
     *  @param array :dataset
     * @return 
     */       
public static void boxNums(Double[] arr) {
        Arrays.sort(arr);
        double Q1 = 0.0;//较小四分位数  Q1的位置=（n+1）/4
        double Q2 = 0.0;//中位数 Q2的位置=（n+1）/2　
        double Q3 = 0.0;//较大四分位数 Q3的位置=3（n+1）/4
        double min = 0.0;//下限  下限=Q1-1.5IQR   四分位距IQR=Q3-Q1   
        double max = 0.0;//上限  上限=Q3+1.5IQR
        List<Double> errNums = new ArrayList<Double>();//异常值   异常值的一个标准：异常值被定义为小于Q1－1.5IQR或大于Q3＋1.5IQR的值 （即大于上限或小于下限的值）
        List<Double> nomalNums = new ArrayList<Double>();//正常值
        Q2 = getMedian(arr);
        if (arr.length % 2 == 0) {//偶数
            Q1 = getMedian(Arrays.copyOfRange(arr, 0, arr.length / 2));
            Q3 = getMedian(Arrays.copyOfRange(arr, arr.length / 2, arr.length));
        } else {//奇数
            Q1 = getMedian(Arrays.copyOfRange(arr, 0, arr.length / 2));
            Q3 = getMedian(Arrays.copyOfRange(arr, arr.length / 2 + 1, arr.length));
        }
        
        System.out.println("q25"+"\t"+Q1);
         System.out.println("q50"+"\t"+Q2);
         System.out.println("q75"+"\t"+Q3);
          System.out.println("max"+"\t"+arr[arr.length-1]);
           System.out.println("min"+"\t"+arr[0]);
     
    }

//求中位数

    public static double getMedian(Double[] arr) {
        if (arr.length == 0) {
            return 0;
        }
        Double[] tempArr = Arrays.copyOf(arr, arr.length);
        Arrays.sort(tempArr);
        if (tempArr.length % 2 == 0) {
            return (tempArr[tempArr.length >> 1] + tempArr[(tempArr.length >> 1) - 1]) / 2.0;
        } else {
            return tempArr[(tempArr.length >> 1)];
        }
    }	

   
        
        
        
}
