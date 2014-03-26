import java.io.File;

/**
 * Created by JSON on 2014/3/24.
 */
public class Test {
    public static void main(String[] arg){

        File folder = new File("/Users/JSON/Documents/MATLAB/num");
        File[] files = folder.listFiles();
        int size = 2701;

        double[][][] x = new double[size-1][][];
        double[][] y = new double[size-1][10];
        for(int i =1; i<size; i++){
            System.out.println("reshape:"+i);
            ImageReader reader = new ImageReader(files[i].getAbsolutePath());
            x[i-1] = reader.reshape(28,28);
            String[] names = files[i].getName().split("-");

            names = names[1].split("\\.");

            int label = Integer.parseInt(names[0])-1;

            y[i-1][label]=1;

        }
        double[][][] train_x =  new double[2000][][];
        double[][] train_y = new double[2000][];
        double[][][] test_x = new double[200][][];
        double[][] test_y = new double[200][];

        System.arraycopy(x,0,train_x,0,2000);
        System.arraycopy(y,0,train_y,0,2000);
        System.arraycopy(x,2500,test_x,0,200);
        System.arraycopy(y,2500,test_y,0,200);

        cnn_layer[] layers = {
                new cnn_layer('i'),
                new cnn_layer('c',"outputmaps","6","kernelsize","5"),
                new cnn_layer('s',"scale","2"),
                new cnn_layer('c',"outputmaps","12","kernelsize","5"),
                new cnn_layer('s',"scale","2")
        };

        cnn_opt opt = new cnn_opt();
        opt.alpha=1;
        opt.batchsize=50;
        opt.numepochs =10;

        cnn_net net = new cnn_net();
        net.layers = layers;
        cnn_net model = cnn.cnnsetup(net, train_x, train_y);

        model = cnn.cnntrain(model, train_x, train_y, opt);
        util.printMatrix(model.o);
        System.out.println("testing");





        System.out.println(cnn.cnntest(model, train_x, train_y));



/*
        double[][] test = {
                {1,1,1,1},
                {1,5,1,1},
                {1,1,1,1},
                {1,1,1,1},
                {0,1,1,1}
        };
        double[][] ker = {
                {1,1},
                {1,1},
                {1,1}
        };
        util.printMatrix(util.convnfull(test,ker));
*/


    }
}
