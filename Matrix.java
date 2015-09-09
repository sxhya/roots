/**
 * Created by user on 1/15/15.
 */
public class Matrix {
    public double[][] myCoo = null;

    public Matrix(double[][] coo) {
        myCoo = new double[coo.length][];
        for (int i = 0; i < myCoo.length; i++) {
            myCoo[i] = new double[coo[i].length];
            System.arraycopy(coo[i], 0, myCoo[i], 0, myCoo[i].length);
        }
    }

    public Matrix(Matrix m) {
        myCoo = new double[m.myCoo.length][];
        for (int i = 0; i < myCoo.length; i++) {
            myCoo[i] = new double[m.myCoo[i].length];
            System.arraycopy(m.myCoo[i], 0, myCoo[i], 0, myCoo[i].length);
        }
    }

    public Matrix(int n) { // Identity matrix constructor
        myCoo = new double[n][n];
        for (int i = 0; i < n; i++) {
            myCoo[i][i] = 1.0;
        }
    }

    int height() {
        return myCoo.length;
    }

    int width() {
        int w = myCoo[0].length;
        for (int i = 1; i < height(); i++)
            if (myCoo[i].length != w)
                throw new IllegalStateException("Not uniform width of matrix rows -- not a matrix");
        return w;
    }

    public Vector mul(Vector v) {
        int w = width();
        if (v.myCoo.length != w) throw new IllegalArgumentException("Matrix width should coincide with vector height");
        double[] result = new double[height()];
        for (int i = 0; i < height(); i++) {
            double c = 0;
            for (int j = 0; j < w; j++) c += myCoo[i][j] * v.myCoo[j];
            result[i] = c;
        }
        return new Vector(result);
    }

    public Matrix mul(Matrix m2) {
        int w = width();
        int w2 = m2.width();
        if (m2.height() != w) throw new IllegalArgumentException("Wrong dimension of matrices");
        double[][] result = new double[height()][w2];
        for (int k = 0; k < w2; k++)
            for (int i = 0; i < height(); i++) {
                double c = 0;
                for (int j = 0; j < w; j++) c += myCoo[i][j] * m2.myCoo[j][k];
                result[i][k] = c;
            }
        return new Matrix(result);
    }

    private static void add(double[][] a, double[] b, int s, int d, double xi) {
        int n = a.length;
        for (int k = 0; k < n; k++) a[d][k] += xi * a[s][k];
        b[d] += xi * b[s];
    }

    public Vector solve(final Vector v) {
        int n;
        if ((n = width()) != height()) throw new IllegalArgumentException();
        Matrix am = new Matrix(this);
        double[][] a = am.myCoo;
        Vector bv = new Vector(v);
        double[] b = bv.myCoo;
        if (n != b.length) throw new IllegalArgumentException();
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            if (a[i][i] == 0) {
                boolean flag = false;
                for (int j = i + 1; j < n; j++)
                    if (a[j][i] != 0) {
                        add(a, b, j, i, 1);
                        flag = true;
                        break;
                    }
                if (!flag) throw new IllegalArgumentException();
            }
            for (int j = i + 1; j < n; j++) add(a, b, i, j, -a[j][i] / a[i][i]);
        }
        for (int i = n - 1; i >= 0; i--) {
            result[i] = b[i];
            for (int j = n - 1; j >= i + 1; j--) result[i] -= result[j] * a[i][j];
            result[i] /= a[i][i];
        }
        return new Vector(result);
    }

    public Matrix transpose() {
        int w;
        int h;
        double[][] coo = new double[w = width()][h = height()];
        for (int i = 0; i < w; i++) for (int j = 0; j < h; j++) coo[i][j] = myCoo[j][i];
        return new Matrix(coo);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Matrix) {
            Matrix m = (Matrix) obj;
            if (m.myCoo.length != myCoo.length) return false;
            for (int i = 0; i < myCoo.length; i++) {
                if (m.myCoo[i].length != myCoo[i].length) return false;
                double norm = 0.0;
                for (int j = 0; j < myCoo[i].length; j++) {
                    double diff = myCoo[i][j] - m.myCoo[i][j];
                    norm += diff * diff;
                }
                if (Math.sqrt(norm) >= Vector.EPS) return false;
            }
            return true;
        }
        return false;
    }

    @Override
    public int hashCode() {
        int result = 0;
        for (int i = 0; i < myCoo.length; i++) {
            for (int j = 0; j < myCoo[i].length; j++) {
                int k = (int) Math.round(myCoo[i][j] / Vector.EPS);
                result = result * 31 + k;
            }
        }
        return result;
    }

    @Override
    public String toString() {
        String result ="";
        for (int i=0; i<myCoo.length; i++) {
            String str = "";
            for (int j = 0; j < myCoo[i].length; j++) {
                String s = Vector.fmt(myCoo[i][j]);
                str += s;
                if (j < myCoo[i].length-1) str+=", ";
            }
            result += str;
            if (i < myCoo.length-1) result+="\n";
        }
        return result;
    }


}
