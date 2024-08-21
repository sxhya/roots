import java.util.*;

class A3Subsystems {
    private fun testA1(a: Vector, b: Vector): Boolean {
        return a == b || a == Vector.minus(b)
    }

    private fun testA2(a: Vector, b: Vector): Boolean {
        if (testA1(a, b)) return true
        val angle = Vector.angle(a, b)
        return Math.abs(angle - 0.5) < Vector.EPS || Math.abs(angle + 0.5) < Vector.EPS
    }

    private fun testA3(a: Vector, b: Vector, rs: Set<Vector>): Int {
        val results = ArrayList<Set<Vector>>()
        if (testA2(a, b)) return -1
        for (v in rs) {
            if (testA2(v, a) && testA2(v, b)) {
                val s = calculateSpan(a, b, v)
                if (s.size != 12) System.err.println("Should be 12 roots")
                var old = false
                for (rs2 in results)
                    if (rs2.containsAll(s)) {
                        old = true
                        break
                    }
                if (!old) results.add(s)
            }
        }
        return results.size
    }

    fun testConjecture(rs: Set<Vector>): Boolean {
        var maxC = 0
        for (a in rs)
            for (b in rs) {
                val a3c = testA3(a, b, rs)
                if (a3c == 0) {
                    println("Found counterexample to A3-conjecture")
                    println("v1 = $a")
                    println("v2 = $b")
                    return false
                } else {
                    if (maxC < a3c) maxC = a3c
                }
            }
        println("No counterexamples found: $maxC maximal No of A3 subsystems")
        return true
    }

    fun testA3Conjecture(a3: Set<Vector>, a: Vector, rootSys: Set<Vector>): Int {
        var sum: MutableSet<Vector> = Utils.unaryOp({ vector -> Vector.add(vector, a) }, a3)

        val intersect = HashSet<Vector>()
        for (r in sum)
            if (rootSys.contains(r)) intersect.add(Vector(r.myCoo))

        if (!intersect.isEmpty()) {
            sum = HashSet(intersect)
            Utils.reflectClosure(sum)
            return if (sum.size == 12) 1 else -1
        }
        return 0
    }

    private fun delta(b: Boolean): Int {
        return if (b) 1 else 0
    }

    private fun findAllA3s(rootSys: Set<Vector>, a: Vector): Set<Set<Vector>> {
        println("Length of a = " + Vector.len(a))
        val a3s = HashSet<Set<Vector>>()
        val c3s = HashSet<Set<Vector>>()
        val b3s = HashSet<Set<Vector>>()
        for (b in rootSys)
            for (c in rootSys)
                if (delta(Vector.orthogonal(a, b)) + delta(Vector.orthogonal(a, c)) + delta(Vector.orthogonal(b, c)) <= 1) {
                    val set = calculateSpan(a, b, c)
                    if (set.size == 12) a3s.add(set)
                    if (set.size == 18 && !c3s.contains(set) && !b3s.contains(set)) {
                        val lC = Utils.lengthCount(set)
                        for (i in lC) print(" $i ")
                        println()
                        if (lC[0] == 6)
                            b3s.add(set)
                        else
                            c3s.add(set)

                    }
                }
        val n: Int
        println("Found #" + a3s.size + " subsystems of type A3")
        println("Found #" + b3s.size + " subsystems of type B3")
        println("Found #" + c3s.size + " subsystems of type C3")
        n = a3s.size + b3s.size + c3s.size

        val a2inc = Array(n) { DoubleArray(n) }
        val r = ArrayList<Set<Vector>>()
        r.addAll(a3s)
        r.addAll(b3s)
        r.addAll(c3s)
        val objects = r.toTypedArray()
        for (i in 0 until n)
            for (j in 0 until n) {
                val rs12 = Utils.intersect(objects[i] as HashSet<Vector>, objects[j] as HashSet<Vector>)
                val size = rs12.size
                when (size) {
                    4, 2, 0, 12, 18 -> a2inc[i][j] = 0.0
                    8 -> {
                        a2inc[i][j] = 2.0
                        a2inc[i][j] = 1.0
                    }
                    6 -> a2inc[i][j] = 1.0
                    else -> throw RuntimeException("Strange intersection : $size")
                }
            }
        for (i in 0 until n) {
            for (j in 0 until n) {
                print(Math.round(a2inc[i][j]).toString() + " ")
            }
            println()
        }

        //Calculate degree of vertices
        /* double deg1 = 0.0;
        for (int i = 0; i < n; i++) {
            double deg = 0;
            for (int j = 0; j < n; j++) {
                deg += a2inc[i][j];
            }
            if (deg1 == 0.0) deg1 = deg;
            else if (deg1 != deg) throw new RuntimeException("Resulting graph not isotropic");
        } */

        var iter = 0
        val a2 = Matrix(a2inc)
        var a2p = Matrix(a2)
        val integers = TreeSet<Long>()
        while (iter++ < 10) {
            a2p = a2p.mul(a2)
            var zeroes = false
            for (i in 0 until n) {
                for (j in 0 until n)
                    if (Math.abs(a2p.myCoo[i][j]) < Vector.EPS) {
                        zeroes = true
                        break
                    } else {
                        integers.add(Math.round(a2p.myCoo[i][j]))
                    }
            }
            if (!zeroes) {
                if (integers.size == 3) {
                    val i = integers.iterator()
                    val mu = i.next()
                    val lambda = i.next()
                    val k = i.next()
                    val m = lambda + 2
                    if (mu == 4L && k == 2 * lambda && (2 * n).toLong() == m * (m - 1))
                        println("Triangular graph T_$m")
                    else
                        println("SRG($n, $k, $lambda, $mu)")
                } else {
                    print("Not a SRG-graph (or there are duplicates in the coefficients): ")
                    for (l in integers) print("$l ")
                    println()
                }
                break
            }
        }
        return a3s
    }

    fun doExperiment(name: String, k: Int, basis: Matrix) {
        var k = k
        val n = basis.height()
        val set = Utils.rootSystemByBasis(basis)
        val tm = basis.transpose()

        var mmr = Vector(basis.myCoo[0])
        val a4s = HashSet<Set<Vector>>()

        println("Total roots (" + name + "): " + set.size)
        k -= 1

        val coos = IntArray(n)
        var min = 0.0
        val lst = ArrayList<Vector>()
        for (v in set) {
            val coo = tm.solve(v)
            var newMin = 0.0
            for (i in 0 until n) newMin += coo.myCoo[i]
            if (newMin < min) {
                mmr = Vector(v.myCoo)
                min = newMin
            }
            for (i in 0 until n)
                if (Math.abs(coo.myCoo[i] - 1) < Vector.EPS) {
                    coos[i]++
                    if (i == k /*&& (v.equals(zero) || set.contains(Vector.minus(v, alphak))) */) {
                        lst.add(v)
                    }
                }
        }
        for (i in 0 until n) print("" + (i + 1) + "-" + coos[i] + ", ")
        println()
        println("Minimal root (wrt to specified order): " + tm.solve(mmr) + " with sum=" + min)
        var pairsCount = 0
        for (i in lst.indices) {
            for (j in i + 1 until lst.size) {
                val rs = HashSet<Vector>()
                rs.add(mmr)
                val vi = lst[i]
                val vj = lst[j]
                rs.add(lst[i])
                rs.add(lst[j])
                Utils.reflectClosure(rs)
                if (rs.size == 6) {
                    //boolean flagi = set.contains(Vector.minus(vi, alphak));
                    //boolean flagj = set.contains(Vector.minus(vj, alphak));
                    println(tm.solve(vi).toString() + "; " + tm.solve(vj))
                    pairsCount++
                }
                if (rs.size == 12) a4s.add(rs)
            }
        }
        println("Constructed " + a4s.size + " subsystems by using unipotent radical (out of " + pairsCount * (pairsCount - 1) / 2 + ")")
        println("Pairs identified: $pairsCount")

        findAllA3s(set, mmr)
        println()
    }

    companion object {
        fun calculateSpan(a: Vector, b: Vector, c: Vector): Set<Vector> {
            val result = HashSet<Vector>()
            result.add(a)
            result.add(b)
            result.add(c)
            Utils.reflectClosure(result)
            return result
        }
    }

}
