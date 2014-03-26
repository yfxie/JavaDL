/**
 * Created by JSON on 2014/3/20.
 */
import java.util.Random;
public class nn {
    public int[] size;
    public int n;
    public double[][][] W;
    public double[][][] vW;
    public double[][][] p;
    public nn(int[] size){
        this.size = size;
        this.n = size.length;
        this.W = new double[this.n-1][][];
        this.vW = new double[this.n-1][][];
        this.p = new double[this.n][1][];
        this.p[0][0] = new double[1];
        for(int i = 1; i<this.n; i++){
            int current_layer_size = this.size[i];
            int prev_layer_size = this.size[i-1];
            this.W[i-1] = new double[current_layer_size][prev_layer_size+1];
            this.vW[i-1] = new double[current_layer_size][prev_layer_size+1];
            this.p[i][0] = new double[current_layer_size];
            // Todo: 其他參數
            // 初始化 Layers
            for(int row = 0 ; row< current_layer_size; row++){
                this.p[i][0][row] = 0.0;
                for(int col=0; col< (prev_layer_size+1); col++){
                    this.W[i-1][row][col] = (Math.random() - 0.5)*2*4*Math.sqrt(6/(current_layer_size+prev_layer_size));
                    this.vW[i-1][row][col] = 0.0;
                }
            }
        }

    }

    public static void main(String[] arg){
        int[] size = {5,4,3};
        nn model = new nn(size);
        util.printMatrix(model.p[2]);

    }

}
