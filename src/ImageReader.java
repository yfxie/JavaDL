import java.awt.*;
import java.awt.image.*;
import java.io.*;
import javax.imageio.*;
public class ImageReader{
    public BufferedImage img;
    public double[] bits;
    public int width;
    public int height;
    public ImageReader(String filename){
        read(filename);
    }
    public double[][] reshape(int x,int y){
        double[][] result = new double[x][y];
        for(int i =0;i<x;i++){
            for(int j=0;j<y;j++){
                result[i][j] = bits[i*x + j%y];
            }
        }
        return  result;
    }
    public void printBits(){
        for(int row = 0; row<height; row++){
            for(int col = 0; col<width; col++){
                System.out.print(bits[width*row+col]+" ");
            }
            System.out.println();
        }
    }
    private void read(String filename){
        try {
            img = ImageIO.read(new File(filename));
            width = img.getWidth();
            height = img.getHeight();
            bits = new double[width*height];
            for(int row = 0; row<height; row++){
                for(int col = 0; col<width; col++){
                    bits[width*row+col] = (img.getRGB(col, row)& 0xFF) / 255.0;
                }
            }
        } catch (IOException e) {
            System.out.println("Error occur.");
        }
    }
}