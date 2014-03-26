import java.util.ArrayList;

/**
 * Created by JSON on 2014/3/21.
 * 注意：Matlab的layer以 w*h*n表示
 *      Java這裡用 n*w*h表示
 *
 */
public class cnn {

    public static cnn_net cnnsetup(cnn_net net, double[][][]train_x, double[][] train_y){
        int inputmaps = 1;

        cnn_layer[] layers = net.layers;
        int one_dim_map_size = train_x[0].length;
        int two_dim_map_size = train_x[0][0].length;
        double[] mapsize = new double[]{one_dim_map_size, two_dim_map_size};
        for(int l = 0 ; l<layers.length; l++){
            if (layers[l].type == 's'){
                mapsize[0] /= layers[l].scale;
                mapsize[1] /= layers[l].scale;
                for (double aMapsize : mapsize) {
                    if (Math.floor(aMapsize) != aMapsize) {
                        System.out.println("Layer " + (l + 1) + " size must be integer. Actual:" + mapsize[0] + " " + mapsize[1]);
                        return null;
                    }
                }
                layers[l].b = new double[inputmaps];
                for(int j = 0; j< inputmaps; j++){
                    layers[l].b[j] = 0;
                }
            }
            if (layers[l].type == 'c'){
                mapsize[0] -= (layers[l].kernelsize-1);
                mapsize[1] -= (layers[l].kernelsize-1);
                int fan_out = layers[l].outputmaps * layers[l].kernelsize * layers[l].kernelsize;
                layers[l].k = new double[inputmaps][layers[l].outputmaps][][];
                layers[l].b = new double[layers[l].outputmaps];
                for(int j = 0; j<layers[l].outputmaps; j++){
                    int fan_in = inputmaps * layers[l].kernelsize * layers[l].kernelsize;
                    for(int i = 0; i<inputmaps;i++){
                        layers[l].k[i][j] = util.multiply(util.diff(util.rand(net.layers[l].kernelsize),0.5),2.0*Math.sqrt(6.0/(fan_in+fan_out)));
                    }
                    layers[l].b[j] = 0;
                }
                inputmaps = layers[l].outputmaps;
            }
        }
        int fvnum = (int)mapsize[0]*(int)mapsize[1]*inputmaps;

        int onum = train_y[0].length;
        net.ffb = new double[onum][1];
        net.ffW = new double[onum][fvnum];
        double ffWFactor = 2.0*Math.sqrt(6.0/(onum+fvnum));

        for(int i = 0; i<onum ;i++){
            for(int j=0; j<fvnum; j++){

                net.ffW[i][j] = (Math.random() - 0.5)*ffWFactor;

            }
        }
        return net;
    }
    public static double cnntest(cnn_net net,double[][][] x,double[][]y){


        net = cnnff(net,x);
        int n = x.length;


        int[] label=new int[n];
        int[] answer = new int[n];
        int bad = 0;

        for(int i =0;i<n;i++){
            double[] tmp = new double[y[0].length];
            for(int k=0;k<tmp.length;k++){
                tmp[k] = net.o[k][i];
            }
            label[i]=util.indexOfMax(tmp);
            answer[i]=util.indexOfMax(y[i]);
            if(label[i]!=answer[i]){
                bad++;
            }
        }
        System.out.println("錯誤數量:"+bad+" 總共:"+n);
        util.printMatrix(net.o);
        return (double)bad/(double)n;
    }
    public static cnn_net cnntrain(cnn_net net,double[][][] train_x,double[][] train_y,cnn_opt opt){
        int m = train_x.length;

        if((m % opt.batchsize) != 0){
            System.out.println("numbatches not integer");
            return null;
        }
        int numbatches = m/opt.batchsize;
        double[][][] x = new double[m][][];
        double[][] y = new double[m][];
        System.arraycopy(train_x,0,x,0,m);
        System.arraycopy(train_y,0,y,0,m);
        y = util.transpose(y);


        net.rL = null;

        for(int i = 0; i < opt.numepochs; i++){
            System.out.println("epoch "+(i+1)+"/"+opt.numepochs);

            int[]kk = util.randperm(m);

            for(int l = 0; l<numbatches; l++){
                int[] partKK = new int[opt.batchsize];

                System.arraycopy(kk,l* opt.batchsize,partKK,0,opt.batchsize);

                double[][][] batch_x = util.copyArray(x,partKK);

                double[][] batch_y = new double[y.length][opt.batchsize];
                for(int y1 =0;y1<y.length;y1++){
                    for(int y2=0;y2<opt.batchsize;y2++){
                        batch_y[y1][y2] = y[y1][partKK[y2]];
                    }
                }

                net = cnnff(net,batch_x);
                net = cnnbp(net,batch_y);
                net = cnnapplygrads(net,opt);
                if(net.rL==null){
                    net.rL =new double[1];
                    net.rL[0] = net.L;
                }
                double[] tmp = net.rL;
                net.rL = new double[net.rL.length+1];
                for(int k=0;k<net.rL.length-1;k++){
                    net.rL[k]=tmp[k];
                }
                net.rL[net.rL.length-1] = 0.99*net.rL[net.rL.length-2]+0.01*net.L;
            }
        }
        return net;
    }
    public static cnn_net cnnapplygrads(cnn_net net,cnn_opt opt){
        for(int l = 1;l<net.layers.length; l++){
            if(net.layers[l].type=='c'){
                for(int j=0;j<net.layers[l].a.length;j++){
                    for(int ii=0;ii<net.layers[l-1].a.length;ii++){
                        net.layers[l].k[ii][j] = util.diff(net.layers[l].k[ii][j], util.multiply(net.layers[l].dk[ii][j],opt.alpha));
                    }
                    net.layers[l].b[j] = net.layers[l].b[j] - net.layers[l].db[j]*opt.alpha;
                }
            }
        }
        net.ffW = util.diff(net.ffW,util.multiply(net.dffW,opt.alpha));
        net.ffb = util.diff(net.ffb,util.multiply(net.dffb,opt.alpha));
        return net;
    }
    public static cnn_net cnnbp(cnn_net net, double[][] y){
        int n =  net.layers.length;

        net.e = util.diff(net.o,y);
        net.L = 1.0/2.0 * util.sum(util.pow(net.e,2.0)) / net.e[0].length;

        net.od = util.dotMutiply(net.e,util.dotMutiply(net.o,util.diff(1.0,net.o)));

        net.fvd = util.multiply(util.transpose(net.ffW),net.od);
        if(net.layers[n-1].type=='c'){
            net.fvd = util.dotMutiply(net.fvd, util.dotMutiply(net.fv,util.diff(1.0,net.fv)) );
        }

        int sa1 = net.layers[n-1].a[0].length;
        int sa2 = net.layers[n-1].a[0][0].length;
        int sa3 = net.layers[n-1].a[0][0][0].length;
        int fvnum = sa2*sa3;
        net.layers[n-1].d = new double[net.layers[n-1].a.length][][][];

        for(int j = 0; j < net.layers[n-1].a.length; j++){
            net.layers[n-1].d[j] = util.reshape( util.subMatrix(net.fvd, j*fvnum ,(j+1)*fvnum,0,net.fvd[0].length), sa1,sa2,sa3);
        }

        for(int l = n-2; l>=0; l--){
            if(net.layers[l].type=='c'){


                net.layers[l].d = new double[net.layers[l].a.length][][][];
                for(int j = 0; j<net.layers[l].a.length; j++){
                    net.layers[l].d[j]= new double[net.layers[l + 1].d[j].length][][];
                    for(int i = 0; i<net.layers[l + 1].d[j].length;i++){
                        double[][] tmp;
                        tmp = util.divide(util.expand(net.layers[l + 1].d[j][i], net.layers[l + 1].scale, net.layers[l + 1].scale),Math.pow(net.layers[l+1].scale,2.0));
                        net.layers[l].d[j][i] = util.dotMutiply( util.dotMutiply(net.layers[l].a[j][i],util.diff(1,net.layers[l].a[j][i])),tmp );
                    }
                }
            }else if(net.layers[l].type=='s'){

                    net.layers[l].d = new double[net.layers[l].a.length][][][];
                    for(int i=0; i <net.layers[l].a.length;i++){
                        double[][][] z = new double[net.layers[l].a[0].length][net.layers[l].a[0][0].length][net.layers[l].a[0][0][0].length];
                        for(int j = 0; j<net.layers[l+1].a.length;j++){
                            z = util.add(z,util.convnfull(net.layers[l + 1].d[j],util.rot180(net.layers[l + 1].k[i][j])));
                        }
                        net.layers[l].d[i]=z;
                    }
            }
        }
        /*  calc gradients */
        for(int l =1; l<n;l++){
            if(net.layers[l].type=='c'){

                net.layers[l].db =new double[net.layers[l].a.length];
                net.layers[l].dk =new double[net.layers[l-1].a.length][net.layers[l].a.length][][];

                for(int j =0;j<net.layers[l].a.length; j++){
                    for(int i =0; i<net.layers[l-1].a.length; i++){

                        net.layers[l].dk[i][j] = util.multiply(util.conv_valid(util.flipall(net.layers[l - 1].a[i]),net.layers[l].d[j]),net.layers[l].d[j].length);
                    }
                    net.layers[l].db[j] = util.sum(net.layers[l].d[j])/net.layers[l].d[j].length;
                }
            }
        }
        net.dffW = util.divide(util.multiply(net.od, util.transpose(net.fv)), net.od[0].length);
        net.dffb = util.mean(net.od,2);

        return net;
    }
    public static cnn_net cnnff(cnn_net net,double[][][] x){

        int n = net.layers.length;
        if(net.layers[0].a==null) {
            net.layers[0].a = new double[1][][][];
        }
        net.layers[0].a[0] = x;
        int inputmaps = 1;
        double[][][] z;
        for(int l = 1; l<n;l++){


            if (net.layers[l].type=='c'){

                net.layers[l].a = new double[net.layers[l].outputmaps][][][];
                for(int j = 0; j<net.layers[l].outputmaps; j++){


                    z  = new double
                            [net.layers[l-1].a[0].length]
                            [net.layers[l-1].a[0][0].length-net.layers[l].kernelsize+1]
                            [net.layers[l-1].a[0][0][0].length-net.layers[l].kernelsize+1];

                    for(int i = 0 ; i<inputmaps; i++){
                        for(int t = 0; t< z.length; t++){
                            z[t] = util.add(z[t],util.conv_valid(net.layers[l - 1].a[i][t], net.layers[l].k[i][j]));
                        }
                    }

                net.layers[l].a[j] = util.sigm(util.add(z,net.layers[l].b[j]));

                }
                inputmaps = net.layers[l].outputmaps;
            }else if(net.layers[l].type == 's'){

                net.layers[l].a = new double[inputmaps][][][];
                for(int j =0;j<inputmaps;j++){

                        z = util.conv_valid(net.layers[l-1].a[j],util.divide(util.ones(net.layers[l].scale), Math.pow(net.layers[l].scale, 2.0)) );
                        int step = net.layers[l].scale;
                        double[][][] tmp = new double[z.length][z[0].length/step+1][z[0][0].length/step+1];
                        for(int i =0;i<z.length;i++){
                            for(int tmp1=0;tmp1<tmp[0].length;tmp1++){
                                for(int tmp2=0;tmp2<tmp[0][0].length;tmp2++){
                                    tmp[i][tmp1][tmp2] = z[i][tmp1*step][tmp2*step];
                                }
                            }
                        }

                        net.layers[l].a[j] = tmp;


                }


            }

        }
        net.fv=null;
        int sa0 = net.layers[n-1].a.length;

        ArrayList<double[][]> tmpList = new ArrayList<double[][]>();
        int fvsize = 0;
        for(int j = 0; j<sa0; j++){
            int sa1 = net.layers[n-1].a[j].length;
            int sa2 = net.layers[n-1].a[j][0].length;
            int sa3 = net.layers[n-1].a[j][0][0].length;

            double[][] tmp = util.reshape(net.layers[n-1].a[j],sa2*sa3,sa1);
            fvsize+=sa2*sa3;
            tmpList.add(tmp);
        }

        net.fv=new double[fvsize][];
        int t=0;
        for (double[][] tmp : tmpList) {
            for (double[] aTmp : tmp) {
                net.fv[t++] = aTmp;
            }
        }

        double[][] net_o_tmp = util.multiply(net.ffW,net.fv);
        double[][] repmat = util.repmat(net.ffb,1,net.fv[0].length);

        net.o = util.sigm(util.add(net_o_tmp,repmat));
        return net;
    }

}
class cnn_net{
    public cnn_layer[] layers;
    public double[][] dffW;
    public double[][] dffb;
    public double[][] ffb;
    public double[][] ffW;
    public double[] rL;
    public double[][] o;
    public double[][] e;
    public double[][] fv;
    public double[][] fvd;
    public double[][] od;
    public double L;
}
class cnn_opt{
    public int batchsize;
    public int numepochs;
    public double alpha;

}
class cnn_layer{
    public double[][][][] dk;
    public double[] db;
    public double a[][][][];
    public double d[][][][];
    public char type;
    public int outputmaps;
    public int kernelsize;
    public int scale;
    public double[] b;
    public double[][][][] k;

    public cnn_layer(char type,String... pairs){
        this.type = type;
        for(int i = 0; i < pairs.length; i+=2){
            int value = Integer.parseInt(pairs[i+1]);
            if(pairs[i].equals("outputmaps")){
                this.outputmaps = value;
            }else if(pairs[i].equals("kernelsize")){
                this.kernelsize = value;
            }else if(pairs[i].equals("scale")){
                this.scale = value;
            }
        }
    }
}