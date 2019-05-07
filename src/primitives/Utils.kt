import java.util.ArrayList
import java.util.HashMap

class Utils {

    companion object {

        fun binaryOp(op: (Pair<Vector, Vector>) -> Vector, a: Set<Vector>, b: Set<Vector>): HashSet<Vector> {
            val result = HashSet<Vector>()
            for (va in a) for (vb in b) result.add(op.invoke(Pair(va, vb)))
            return result
        }

        fun unaryOp(op: (Vector) -> Vector, a: Set<Vector>): HashSet<Vector> {
            val result = HashSet<Vector>()
            for (va in a) result.add(op.invoke(va))
            return result
        }

        private fun reflect(set: Set<Vector>, wrt: Set<Vector>) = binaryOp({ p: Pair<Vector, Vector> -> p.first.reflectWRT(p.second) }, wrt, set)

        private fun reflectJoin(set: MutableSet<Vector>, wrt: Set<Vector>) {
            set.addAll(reflect(set, wrt))
        }

        fun reflectClosure(set: MutableSet<Vector>) {
            var oldLength: Int
            do {
                oldLength = set.size
                reflectJoin(set, set)
            } while (set.size > oldLength)
        }

        fun rootSystemByBasis(basis: Matrix): HashSet<Vector> {
            val result = HashSet<Vector>()
            for (i in 0 until basis.myCoo.size) result.add(Vector(basis.myCoo[i]))
            reflectClosure(result)
            return result
        }

        fun<T> intersect(a: Set<T>, b: Set<T>): Set<T> = a.intersect(b)

        fun lengthCount(myRoots: Set<Vector>): List<Int> {
            val lengths = ArrayList<Double>()
            val counter = HashMap<Double, Int>()
            for (v in myRoots) {
                val vecLen = Vector.len(v)
                var key: Double = vecLen
                for (k in lengths) if (Math.abs(k - vecLen) < Vector.EPS) key = k
                var s: Int? = counter[key]
                if (s == null) lengths.add(key)
                s = if (s == null) 1 else s + 1
                counter[key] = s
            }
            lengths.sort()
            val result = ArrayList<Int>()
            for (i in lengths) result.add(counter[i]!!)

            return result
        }

        fun subtractClosure(set: Set<Vector>, brs: BasedRootSystem): HashSet<Vector> {
            val result = HashSet(set)
            var counter = 0
            do {
                val extraRoots = HashSet<Vector>()
                for (v1 in result) for (v2 in result) {
                    val element = Vector.minus(v1, v2)
                    if (brs.myRootSet.contains(element))
                        if (!result.contains(element))
                            extraRoots.add(element)
                }
                val srs = brs.myFundamentalRoots.filter { result.contains(it) }.size
                val srm = brs.myFundamentalRoots.filter { result.contains(Vector.minus(it)) }.size
                System.out.println("Iteration: $counter Roots: ${result.size} SimpleRoots: $srs MinusSimpleRoots: $srm")
                result.addAll(extraRoots)
                counter++
            } while (extraRoots.isNotEmpty())
            System.out.println("Counter: $counter")
            return result
        }

    }

    }