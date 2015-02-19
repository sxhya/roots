import java.util.*;
import java.util.concurrent.Callable;

/**
 * Created by sxh on 11/8/14.
 */
public class RootSet {

    public abstract static class VectorUOp {
        public abstract Vector exec (Vector v);
    }

    public abstract static class VectorBOp {
        public abstract Vector exec (Vector a, Vector b);
    }


    public Set<Vector> myRoots = new HashSet<Vector>();

    public RootSet(){}

    public RootSet (Set roots) {
        myRoots = roots;
    }

    public RootSet(Matrix basis) {
        for (int i=0; i< basis.myCoo.length; i++) addRoot(new Vector(basis.myCoo[i]));
        reflectClosure();
    }

    public RootSet(RootSet set) {
        for (Vector v : set.myRoots)
          myRoots.add(new Vector(v));
    }

    public void addRoot(Vector r) {
        boolean flag = false;
        for (Vector r2 : myRoots) if (r2.equals(r)) {flag = true; break;}
        if (!flag) myRoots.add(r);
    }

    public boolean contains(Vector r) {
        return myRoots.contains(r);
    }

    public boolean contains(RootSet rs) {
        for (Vector v : rs.myRoots) if (!contains(v)) return false;
        return true;
    }

    public static RootSet join (final RootSet rs1, final RootSet rs2) {
        RootSet result = new RootSet(rs1); result.join(rs2);
        return result;
    }

    public void join (final RootSet rs1) {
        RootSet rs2 = new RootSet(rs1);
        for (Vector v : rs2.myRoots) addRoot(v);
    }

    public static RootSet binaryOp(VectorBOp op, final RootSet a, final RootSet b) {
        RootSet result = new RootSet();
        for (Vector va : a.myRoots) for (Vector vb : b.myRoots) {
            result.addRoot(op.exec(va, vb));
        }
        return result;
    }

    public static RootSet unaryOp(VectorUOp op, RootSet a) {
        RootSet result = new RootSet();
        for (Vector va : a.myRoots) {
            result.addRoot(op.exec(va));
        }
        return result;
    }

    public static RootSet reflectJoin(RootSet set, RootSet wrt) {
        RootSet result = new RootSet(set);
        result.reflectJoin(wrt);
        return result;
    }

    public void reflectJoin(RootSet wrt) {
        join (binaryOp(new VectorBOp() {
            @Override
            public Vector exec(Vector a, Vector b) {
                return a.reflectWRT(b);
            }
        }, wrt, this));
    }


    public static RootSet mirrorJoin(RootSet set) {
        RootSet result = new RootSet(set);
        RootSet reflected = unaryOp(new VectorUOp() {
            @Override
            public Vector exec(Vector a) {
                return Vector.minus(a);
            }
        }, set);
        return join (result, reflected);
    }

    public RootSet reflectClosure () {
        int oldLength;
        do {oldLength = myRoots.size(); reflectJoin(this);
        } while (myRoots.size() > oldLength);
        return this;
    }

    public RootSet perpendicularTo(Vector v) {
        RootSet result = new RootSet();
        for (Vector v2 : myRoots) {
            if (Math.abs(Vector.prod(v ,v2))<Vector.EPS) result.myRoots.add(v2);
        }
        return result;
    }

    public static Set<Vector> intersect (Set<Vector> s1, Set<Vector> s2) {
        Set<Vector> result = new HashSet<Vector>();
        for (Vector v : s1)
         if (s2.contains(v)) result.add(v);
        return result;
    }
}
