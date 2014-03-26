import java.util.Random;

/**
 * Created by JSON on 2014/3/20.
 */
public class util {
    public static void main(String[] arg){
        double[][] test={
                {0,1},{2,3}


        };
        double[][][] a= reshape(test,1,2,2);
        printMatrix(a[0]);
    }
    static double[][] divide(double[][] a,double b){
        int row =a.length;
        double[][] result = new double[row][];
        for(int i=0;i<row;i++)
            result[i]=divide(a[i],b);
        return result;
    }
    static double[][][] divide(double[][][]a,double b){
        int row=a.length;
        double[][][] result = new double[row][][];
        for(int i=0;i<row;i++)
            result[i]=divide(a[i],b);
        return result;
    }
    static double[] divide(double[] a,double b){
        int row = a.length;
        double[] result = new double[row];
        for(int i =0;i<row;i++)
            result[i] = a[i]/b;
        return result;
    }
    public static int indexOfMax(double[]a){
        int max = 0;
        for(int i =1;i<a.length;i++){
            if(a[i]>a[max])
                max = i;
        }
        return max;
    }
    public static  void printMatrix(int[] m){
        for(int i =0;i<m.length;i++)
            System.out.print(m[i] + " ");

        System.out.println();
    }
    public static void printMatrix(double[][] m){
        for(int i = 0; i < m.length; i++){
            for(int k=0; k < m[i].length; k++){
                System.out.print(m[i][k]+" ");
            }
            System.out.println();

        }
    }
    public static int[] randperm(int m){
        int[] array = new int[m];
        for(int i = 0; i< m ; i++){
            array[i] = i;
        }
        shuffleArray(array);
        return array;
    }
    static double[][] rot180(double[][] a){
        return flipdim(flipdim(a,1),2);
    }
    static double[][] flipdim(double[][] a,int deg){
        int row = a.length;
        int col = a[0].length;
        double[][] result = new double[row][col];
        if(deg==1){
            for(int i=0;i<row;i++){
                result[i] = a[row-i-1];
            }
        }else if(deg==2){
            for(int i=0;i<row;i++){
                for(int j =0;j<col;j++){
                    result[i][j] = a[i][col-j-1];
                }
            }
        }
        return result;
    }
    static double[][] reshape(double[][][] a,int x, int y){
        int a1 = a.length;
        int a2 = a[0].length;
        int a3 = a[0][0].length;
        int size = x*y;

        if(a1*a2*a3 != size){
            System.out.println("大小不一致");
            return null;
        }
        double[][] output = new double[x][y];
        for(int i =0; i<size; i++){

            output[i%x][i/x] = a[i/(a2*a3)][(i/a3)%a2][i%a3];
        }
        return output;
    }
    // Ex: repmat([1 2;3 4],1,2) => [1 2 1 2;3 4 3 4]
    static double[][] repmat(double[][] a,int x,int y){
        double[][] result = new double[a.length*x][a[0].length*y];
        for(int i = 0; i<result.length; i++){
            for(int j=0;j<result[0].length;j++){
                int row = i % a.length;
                int col = j % a[0].length;
                result[i][j]= a[row][col];
            }
        }
        return result;
    }
    static double[][] dotMutiply(double a[][],double[][] b){

        int row = a.length;
        int col = a[0].length;
        int rowb=b.length;
        int colb=b[0].length;
        if ( row != rowb || col!=colb ) {
            throw new IllegalArgumentException("大小不一致");
        }
        double[][] result = new double[row][col];
        for(int i =0; i<row; i++){
            for(int j=0;j<col;j++){
                result[i][j] = a[i][j]*b[i][j];
            }
        }
        return result;
    }
    static  double[][][] multiply(double[][][] a,double b){
        int row = a.length;
        double[][][] result = new double[row][][];
        for(int i =0;i<row;i++){
            result[i] = multiply(a[i],b);

        }
        return result;
    }
    static double[] multiply(double[] a, double b){
        int col = a.length;
        double[] result = new double[col];
        for(int i = 0;i<col;i++){
            result[i] = a[i]*b;
        }
        return  result;
    }
    static double[][] multiply(double a[][],double b){
        int row = a.length;
        int col = a[0].length;
        double[][] result = new double[row][];
        for(int i =0; i<row; i++){
            result[i] = multiply(a[i],b);
        }
        return result;
    }
    static double[][] multiply(double a[][], double b[][]) {

        int aRows = a.length,
                aColumns = a[0].length,
                bRows = b.length,
                bColumns = b[0].length;

        if ( aColumns != bRows ) {
            throw new IllegalArgumentException("A:Rows: " + aColumns + " did not match B:Columns " + bRows + ".");
        }

        double[][] resultant = new double[aRows][bColumns];

        for(int i = 0; i < aRows; i++) { // aRow
            for(int j = 0; j < bColumns; j++) { // bColumn
                for(int k = 0; k < aColumns; k++) { // aColumn
                    resultant[i][j] += a[i][k] * b[k][j];
                }
            }
        }

        return resultant;
    }

    static double[][] diff(double a,double[][] b){
        int row = b.length;
        int col = b[0].length;
        double[][] result = new double[row][col];
        for(int i =0; i<row; i++){
            for(int j=0;j<col;j++){
                result[i][j] = a - b[i][j];
            }
        }
        return result;
    }
    static double[][] expand(double[][] a, int x, int y){
        double[][] result = new double[a.length*x][a[0].length*y];
        for(int i = 0; i < result.length; i++){
            for(int j =0; j<result[0].length;j++){
                result[i][j] = a[i/x][j/y];
            }
        }
        return result;
    }

    static double[][] subMatrix(double[][] a,int row_start,int row_end, int col_start,int col_end){
        //System.out.println(row_start+" "+row_end+" "+col_start+" "+col_end);
        double[][] result = new double[row_end-row_start][col_end-col_start];
        for(int i=row_start; i <row_end; i++){
            for(int j = col_start; j<col_end;j++){
                result[i-row_start][j-col_start] = a[i][j];
            }
        }
        return result;
    }

    static double[][][] reshape(double[][] a,int level1, int level2, int level3){

        double[][][] result = new double[level1][level2][level3];
        int row = a.length;
        int col = a[0].length;
        if ( row*col != level1*level2*level3 ) {
            throw new IllegalArgumentException("大小不一致.");
        }
        int size = row*col;
        for(int i =0;i<size;i++){
            result[i/(level2*level3)][(i/level3)%level2][i%level3]=a[i%row][i/row];
        }
        return result;
    }

    static double[][] transpose(double[][] a){
        double[][] result = new double[a[0].length][a.length];
        for(int i =0 ; i<a.length;i++){
            for(int j = 0;j<a[0].length;j++){
                result[j][i] = a[i][j];
            }
        }
        return result;
    }
    static double[] diff(double[] a, double[] b){
        int row = a.length;
        double[] result = new double[row];
        for(int i =0;i<row;i++){
            result[i] = a[i]-b[i];
        }
        return result;
    }
    static double[][] diff(double[][] a,double b){
        int row = a.length;
        int col = a[0].length;
        double[][] result = new double[row][col];
        for(int i =0; i<row; i++){
            for(int k=0;k<col;k++){
                result[i][k] = a[i][k]-b;
            }
        }
        return result;
    }
    static double[][] diff(double[][] a, double[][] b){
        int row = a.length;
        int col = a[0].length;
        int rowb = b.length;
        int colb = b[0].length;
        if(row!=rowb || col!=colb){
            System.out.println("大小不一致");
            return null;
        }
        double[][] result = new double[row][col];
        for(int i =0; i<row; i++){
            result[i] = diff(a[i],b[i]);
        }
        return result;
    }
    static double[][] rand(int size){
        double[][] result = new double[size][size];
        for(int i=0;i<size;i++){
            for(int k=0;k<size;k++)
                result[i][k] = Math.random();
        }
        return result;
    }
    static double[][] mean(double[][] a,int type){
        double[][] result;
        int row=a.length;
        int col=a[0].length;
        if(type==1){
            result = new double[1][col];
            for(int i =0 ; i<col; i++){
                for(int j =0; j<row; j++){
                    result[0][i] += a[j][i];
                }
                result[0][i] /= row;
            }
            return result;
        }else if(type==2){
            result = new double[row][1];
            for(int i=0; i<row; i++){
                for(int j=0;j<col;j++){
                    result[i][0] += a[i][j];
                }
                result[i][0]/=col;
            }
            return result;
        }
        return null;
    }

    static double sum(double[][][] a){
        double result = 0.0;
        int row = a.length;
        for(int i=0; i < row; i++){
            result+=sum(a[i]);
        }
        return result;
    }
    static double sum(double[][] a){
        int row = a.length;
        int col = a[0].length;
        double result = 0.0;
        for(int i = 0; i<row; i++){
            for(int j=0; j<col; j++){
                result += a[i][j];
            }
        }
        return result;
    }
    static double[][] pow(double[][] a, double p){
        int row = a.length;
        int col = a[0].length;
        double[][] result = new double[row][col];
        for(int i = 0; i<row; i++){
            for(int j=0; j<col; j++){
                result[i][j] = Math.pow(a[i][j],p);
            }
        }
        return result;
    }
    static double sigm(double p){
        return 1.0/(1.0+Math.exp(-p));
    }
    static double[][] sigm(double[][] a){

        int row =a.length;
        int col = a[0].length;
        double[][] result = new double[row][col];
        for(int i=0;i<row;i++){
            for(int j =0;j<col;j++){
                result[i][j] = sigm(a[i][j]);
            }
        }
        return result;
    }
    static double[][][] sigm(double[][][] a){
        int a1 = a.length;
        int a2 = a[0].length;
        int a3 = a[0][0].length;
        double[][][] result = new double[a1][a2][a3];
        for(int i =0; i<a1;i++){
            for(int j=0;j<a2;j++){
                for(int k=0;k<a3;k++){
                    result[i][j][k] = sigm(a[i][j][k]);
                }
            }
        }
        return result;
    }
    static double[][][] copyArray(double[][][] x,int[] index){
        double[][][] result = new double[index.length][][];
        for(int i =0;i<result.length;i++){
            result[i] = x[index[i]];
        }
        return result;
    }
    static double[][] copyArray(double[][] x,int[] index){
        double[][] result = new double[index.length][];
        for(int i =0;i<result.length;i++){
            result[i] = x[index[i]];
        }
        return result;
    }
    static void shuffleArray(int[] ar)
    {
        Random rnd = new Random();
        for (int i = ar.length - 1; i > 0; i--)
        {
            int index = rnd.nextInt(i + 1);
            // Simple swap
            int a = ar[index];
            ar[index] = ar[i];
            ar[i] = a;
        }
    }
    public static double singlePixelConvolution(double[][] input,
                                                int x, int y,
                                                double[][] k,
                                                int kernelWidth,
                                                int kernelHeight) {
        double output = 0;
        for (int i = 0; i < kernelWidth; ++i) {
            for (int j = 0; j < kernelHeight; ++j) {
                output = output + (input[x + i][y + j] * k[i][j]);
            }
        }
        return output;
    }
    static void printForMatlab(double[][][] train_x){
        int level1= train_x.length;
        int level2= train_x[0].length;
        int level3= train_x[0][0].length;
        for(int i=0;i<level1*level2*level3;i++){
            System.out.print(train_x[i / (level2 * level3)][(i / level3) % level2][i % level3] + " ");
        }
    }
    static double[][] flipall(double[][] a){
        return flipdim(flipdim(a, 1), 2);
    }
    static double[][][] flipall(double[][][] a){

        double[][][] result = new double[a.length][][];
        for(int i=0;i<a.length;i++){
            result[a.length-i-1] = flipall(a[i]);
        }
        return result;
    }
    static double[][][] convnfull(double[][][] a,double[][] kernel){
        double[][][] result = new double[a.length][][];
        for(int i =0;i<result.length;i++)
            result[i] = convnfull(a[i],kernel);
        return result;
    }
    static double[][] convnfull(double[][] a, double[][] kernel){
        double[][] result  = new double[a.length+(kernel.length-1)*2][a[0].length+(kernel[0].length-1)*2];
        double[][] result2  = new double[a.length+(kernel.length-1)][a[0].length+(kernel[0].length-1)];
        for(int i =0;i<a.length;i++){
            for(int j=0;j<a[i].length;j++){

                result[kernel.length-1+i][kernel[0].length-1+j] =a[i][j];
            }
        }
        for(int i =0;i<result.length;i++){
            for(int j =0;j<result[i].length;j++){
                double sum=0.0;
                for(int k1=0;k1<kernel.length;k1++){
                    for(int k2=0;k2<kernel[0].length;k2++){
                        int x = k1+i;
                        int y = k2+j;
                        if(x<result.length && y<result[0].length){
                            sum += result[x][y]*kernel[k1][k2];
                        }
                    }
                }
                if(i<result2.length && j<result2[0].length)
                    result2[i][j]=sum;
            }
        }


        return result2;
    }
    static double[][] add(double[][] a,double[][] b){
        int a1 = a.length;
        int a2 = a[0].length;
        int b1 = b.length;
        int b2 = b[0].length;
        if(a1!=b1 || a2!=b2){
            throw new IllegalArgumentException("a:"+a1+"x"+a2+" b:"+b1+"x"+b2);

        }
        double[][] result = new double[a1][a2];
        for(int i=0;i<a1;i++){
            for(int j=0;j<a2;j++){
                result[i][j] = a[i][j]+b[i][j];
            }
        }
        return result;
    }
    static double[][][] add(double[][][]a,double b){
        int a1 = a.length;
        int a2 = a[0].length;
        int a3 = a[0][0].length;
        double[][][] result = new double[a1][a2][a3];
        for(int i =0; i<a1;i++){
            for(int j=0;j<a2;j++){
                for(int k=0;k<a3;k++){
                    result[i][j][k] = a[i][j][k]+b;
                }
            }
        }
        return result;
    }
    static double[][][] add(double[][][]a,double[][][]b){
        int a1 = a.length;
        int a2 = a[0].length;
        int a3 = a[0][0].length;
        double[][][] result = new double[a1][a2][a3];
        for(int i =0; i<a1;i++){
            for(int j=0;j<a2;j++){
                for(int k=0;k<a3;k++){
                    result[i][j][k] = a[i][j][k]+b[i][j][k];
                }
            }
        }
        return result;
    }
    static double[][][] conv_valid(double[][][] input, double[][] kernel){
        double[][][] result = new double[input.length][][];
        for(int i =0; i<input.length;i++){
            result[i] = conv_valid(input[i], kernel);
        }
        return result;
    }

    static double[][] conv_valid(double[][][] input,double[][][] kernel){
        double[][] result = null;
        for(int i =0;i<input.length;i++){
            if(result==null){

                result = conv_valid(input[i], kernel[i]);
            }else{

                result = add(result,conv_valid(input[i], kernel[i]));
            }
        }
        return result;
    }
    static double[][] conv_valid(double[][] input,double[][] kernel) {
        int width = input.length;
        int height = input[0].length;
        int kernelWidth = kernel[0].length;
        int kernelHeight = kernel.length;
        int smallWidth = width - kernelWidth + 1;
        int smallHeight = height - kernelHeight + 1;
        double[][] output = new double[smallWidth][smallHeight];
        for (int i = 0; i < smallWidth; ++i) {
            for (int j = 0; j < smallHeight; ++j) {
                output[i][j] = 0;
            }
        }
        for (int i = 0; i < smallWidth; ++i) {
            for (int j = 0; j < smallHeight; ++j) {
                output[i][j] = singlePixelConvolution(input, i, j, kernel,
                        kernelWidth, kernelHeight);
            }
        }
        return output;
    }

    static double[][] ones(int size){
        double[][] result = new double[size][size];
        for(int i = 0; i<size ;i++){
            for(int j=0; j<size; j++){
                result[i][j] = 1.0;
            }
        }
        return result;
    }

    static void shuffleArray(double[][] ar)
    {
        Random rnd = new Random();
        for (int i = ar.length - 1; i > 0; i--)
        {
            int index = rnd.nextInt(i + 1);
            // Simple swap
            double[] a = ar[index];
            ar[index] = ar[i];
            ar[i] = a;
        }
    }
    static void shuffleArray(double[][][] ar)
    {
        Random rnd = new Random();
        for (int i = ar.length - 1; i > 0; i--)
        {
            int index = rnd.nextInt(i + 1);
            // Simple swap
            double[][] a = ar[index];
            ar[index] = ar[i];
            ar[i] = a;
        }
    }
}
