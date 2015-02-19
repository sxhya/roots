/**
 * Created by sxh on 11/8/14.
 */
public class Vector{
    public static final double EPS = 0.01;

    public double[] myCoo = null;

    public Vector () {}

    public Vector (double[] coo) {
        myCoo = coo;
    }

    public Vector (Vector c) {
        myCoo = new double[c.myCoo.length];
        for (int i=0; i<myCoo.length; i++) myCoo[i] = c.myCoo[i];
    }

    public static Vector add (Vector a, Vector b) {
        if (a.myCoo.length != b.myCoo.length) throw new IllegalArgumentException();
        Vector result = new Vector();
        result.myCoo = new double[a.myCoo.length];
        for (int i=0; i< a.myCoo.length; i++) result.myCoo[i] = a.myCoo[i] + b.myCoo[i];
        return result;
    }

    public static Vector mul (Vector a, double c) {
        Vector result = new Vector();
        result.myCoo = new double[a.myCoo.length];
        for (int i=0; i< a.myCoo.length; i++) result.myCoo[i] = c * a.myCoo[i];
        return result;
    }

    public static Vector minus (Vector a) {return mul(a, -1);}

    public static Vector minus (Vector a, Vector b) {return add(a, minus(b));}

    public static double len2 (Vector a) {
        double result = 0.0;
        for (int i=0; i< a.myCoo.length; i++) result += a.myCoo[i]*a.myCoo[i];
        return result;
    }

    public static double len (Vector a) {return Math.sqrt(len2(a));}

    public static double prod (Vector a, Vector b) {
        if (a.myCoo.length != b.myCoo.length) throw new IllegalArgumentException();
        double result = 0.0;
        for (int i=0; i< a.myCoo.length; i++) result += a.myCoo[i]*b.myCoo[i];
        return result; /// (len(a) * len(b));
    }

    public static double angle (Vector a, Vector b) {
        return Vector.prod(a, b) / (Vector.len(a) * Vector.len(b));
    }

    public static boolean orthogonal(Vector a, Vector b) {
        return Math.abs(angle(a, b)) < EPS;
    }

    public Vector reflectWRT (Vector y) {return add(y, mul(this, - 2 * prod (this, y) / len2(this)));}

    @Override
    public boolean equals(Object o) {
        if (o == null) return false;
        if (o instanceof Vector) {
            Vector vo = (Vector) o;
            return (len(add(this, minus(vo))) < EPS);
        }
        throw new IllegalArgumentException();
    }

    @Override
    public int hashCode() {
        int result = 0;
        for (double aMyCoo : myCoo) {
            int k = (int) Math.round(aMyCoo / EPS);
            result = result * 31 + k;
        }
        return result;
    }

    @Override
    public String toString() {
        String result = "{";
        for (int i=0; i<myCoo.length; i++) {
            result += fmt(myCoo[i]);
            if (i < myCoo.length-1) result+=", ";
        }
        return result + "}";
    }

    public static String fmt(double d) {
        if(Math.abs((d - (long) d)) <= 0.00001) return String.format("%d",(long)d); else return String.format("%e",d);}
}
