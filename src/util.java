import java.util.Random;

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
    public static void printMatrix(double[] m){
        for(int i=0;i<m.length;i++)
            System.out.print(m[i]+" ");
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
    static double[][] rot180(double[][] a){ //正確
        return flipdim(flipdim(a,1),2);
    }
    static double[][] flipdim(double[][] a,int deg){ //正確
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

            output[i%x][i/x] = a[i/(a2*a3)][i%a2][(i/a2)%a3];
        }
        return output;
    }
    // Ex: repmat([1 2;3 4],1,2) => [1 2 1 2;3 4 3 4]
    static double[][] repmat(double[][] a,int x,int y){
        double[][] result = new double[a.length*x][a[0].length*y];
        for(int i = 0; i<result.length; i++){
            for(int j=0;j<result[i].length;j++){
                int row = i % a.length;
                int col = j % a[0].length;
                result[i][j]= a[row][col];
            }
        }
        return result;
    }
    static double[][][] dotMutiply(double a[][][], double[][][] b){
        int rowa = a.length;
        int rowb = b.length;
        if(rowa != rowb){
            throw new IllegalArgumentException("大小不一致");
        }
        double[][][] result = new double[rowa][][];
        for(int i =0;i<rowa;i++){
            result[i] = dotMutiply(a[i],b[i]);
        }
        return result;
    }
    static double[][] dotMutiply(double a[][],double[][] b){
        int rowa = a.length;
        int rowb = b.length;
        if(rowa != rowb){
            throw new IllegalArgumentException("大小不一致");
        }
        double[][] result = new double[rowa][];
        for(int i =0;i<rowa;i++){
            result[i] = dotMutiply(a[i],b[i]);
        }
        return result;
    }
    static  double[] dotMutiply(double[] a,double[] b){
        int rowa = a.length;
        int rowb = b.length;
        if(rowa != rowb){
            throw new IllegalArgumentException("大小不一致");
        }
        double[] result = new double[rowa];
        for(int i =0;i<rowa;i++){
            result[i] = a[i]*b[i];
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
    static double[][][] diff(double a,double[][][] b){
        int row = b.length;
        double[][][] result = new double[row][][];
        for(int i =0; i<row; i++){
            result[i] = diff(a,b[i]);
        }
        return result;

    }
    static double[][] diff(double a,double[][] b){
        int row = b.length;
        double[][] result = new double[row][];
        for(int i =0; i<row; i++){
            result[i] = diff(a,b[i]);
        }
        return result;
    }
    static double[] diff(double a,double[] b){
        int row = b.length;
        double[] result = new double[row];
        for(int i=0;i<row;i++){
            result[i] = a-b[i];
        }
        return result;
    }
    static double[][][] expand(double[][][] a,int x,int y){
        int row = a.length;
        double[][][] result = new double[row][][];
        for(int i =0;i<row;i++){
            result[i] = expand(a[i],x,y);
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
            result[i/(level2*level3)][(i/level3)%level2][i%level3] = a[i/col][i%col];
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
        double result = 0.0;
        int row = a.length;
        for(int i = 0; i<row; i++){
            result += sum(a[i]);
        }
        return result;
    }
    static double sum(double[] a){
        double result = 0.0;
        int row=a.length;
        for(int i=0;i<row;i++){
            result+=a[i];
        }
        return result;
    }
    static double[][] dotPow(double[][] a, double p){
        int row = a.length;
        double[][] result = new double[row][];
        for(int i = 0; i<row; i++){
            result[i] = new double[a[i].length];
            for(int j=0; j<a[i].length; j++){
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
    static double[][][] convnfull(double[][][] a,double[][] kernel){ //正確
        double[][][] result = new double[a.length][][];
        for(int i =0;i<result.length;i++)
            result[i] = convnfull(a[i],kernel);
        return result;
    }
    static double[][] convnfull(double[][]a,double[][] kernel){
        int krow = kernel.length;
        int arow = kernel.length;
        double[][] result = new double[a.length+kernel.length-1][a[0].length+kernel[0].length-1];
        double[][][] tmp_result = new double[kernel.length][][];
        for(int i =0;i<tmp_result.length;i++){
            tmp_result[i] = convnfull(a,kernel[i]);
        }
        for(int i =0;i<tmp_result.length;i++){
            for(int j=0;j<tmp_result[i].length;j++){
                for(int t=0;t<tmp_result[i][j].length;t++){
                    result[i+j][t]+=tmp_result[i][j][t];
                }
            }
        }
        return result;
    }
    static double[][] convnfull(double[][] a,double[] kernel){
        double[][] result = new double[a.length][];
        for(int i =0;i<result.length;i++){
            result[i] = convnfull(a[i],kernel);
        }
        return result;
    }
    static double[] convnfull(double[] a,double[] kernel){

        double[] result = new double[a.length+kernel.length-1];
        for(int i=0;i<kernel.length;i++){
            for(int j =0;j<a.length;j++){
                result[j+i] += a[j]*kernel[i];
            }
        }
        return result;
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
    static double[][] conv_valid(double[][][] a,double[][][] b){

        double[][][] result = new double[a.length-b.length+1][a[0].length-b[0].length+1][a[0][0].length-b[0][0].length+1];
        for(int r =0;r<result.length;r++){
            for(int c=0;c<result[0].length;c++){
                for(int z=0;z<result[0][0].length;z++){
                    for(int rb=b.length-1;rb>=0;rb--){
                        for(int cb=b[0].length-1;cb>=0;cb--){
                            for(int cz=b[0][0].length-1;cz>=0;cz--)
                                result[r][c][z] = result[r][c][z]+a[r+b.length-rb-1][c+b[0].length-cb-1][z+b[0][0].length-cz-1]*b[rb][cb][cz];
                        }
                    }
                }

            }
        }

        return result[0];

    }
    static double[][][] conv_valid(double[][][] a,double[][] b){
        double[][][] result = new double[a.length][][];
        for(int i=0;i<a.length;i++)
            result[i]=conv_valid(a[i],b);
        return result;
    }
    static double[][] conv_valid(double[][] a,double[][] b){
        double[][] result = new double[a.length-b.length+1][a[0].length-b[0].length+1];
        for(int r =0;r<result.length;r++){
            for(int c=0;c<result[0].length;c++){
                for(int rb=b.length-1;rb>=0;rb--){
                    for(int cb=b[0].length-1;cb>=0;cb--){
                        result[r][c] = result[r][c]+a[r+b.length-rb-1][c+b[0].length-cb-1]*b[rb][cb];
                    }
                }
            }
        }
        return result;
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
