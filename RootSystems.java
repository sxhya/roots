/**
 * Created by user on 1/17/15.
 */
public class RootSystems {

    public static final Matrix e8base = new Matrix(new double[][]{
            {-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5},
            { 0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  0.0},
            { 0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  0.0},
            { 0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0},
            { 0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  0.0},
            { 0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0},
            { 0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0},
            { 1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}});

    public static final Matrix e7base = new Matrix(new double[][]{
            {-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, Math.sqrt(2)/2},
            {   0,    0,    0,    0,    1,   -1, 0},
            {   0,    0,    0,    0,    1,    1, 0},
            {   0,    0,    0,    1,   -1,    0, 0},
            {   0,    0,    1,   -1,    0,    0, 0},
            {   0,    1,   -1,    0,    0,    0, 0},
            {   1,   -1,    0,    0,    0,    0, 0}});

    public static final Matrix e6base = new Matrix(new double[][]{
            {   1,   -1,   0,     0,    0, 0},
            {   0,    0,    0,    1,   -1, 0},
            {   0,    1,   -1,    0,    0, 0},
            {   0,    0,    1,   -1,    0, 0},
            {   0,    0,    0,    1,    1, 0},
            {-0.5, -0.5, -0.5, -0.5, -0.5, Math.sqrt(3)/2}});

    public static Matrix dlbase (int l) {
        double[][] coos = new double[l][l];
        for (int i=0; i<l; i++)  {
            coos[i][i] = 1.0;
            if (i < l-1) coos[i][i+1] = -1.0;}
        coos[l-1][l-2] = 1.0;
        return new Matrix(coos);
    }

    public static Matrix albase (int l) {
        double[][] coos = new double[l][l];
        for (int i=0; i<l-1; i++)  {
            coos[i][i] = 1.0;
            coos[i][i+1] = -1.0;}
        double b = ((l-1) + Math.sqrt(l+1))/l;
        for (int i=0; i<l-1; i++) coos[l-1][i] = b - 1;
        coos[l-1][l-1] = b;
        return new Matrix(coos);
    }
}
