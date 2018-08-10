import java.util.*

/**
 * Created by sxh on 8/28/15.
 */

//TODO: I don'like the arguments of the primary constructor...
class WeylGroup (basis: Matrix, cutoff: Int, noReflection: Int) {//if you want to get full Weyl group and not the Levi Weyl subgroup noReflection should be -1
    val carrier: MutableSet<Matrix> = HashSet()

    companion object {
        fun lol(brs: BasedRootSystem, sim: Set<Vector>): Pair<Matrix, Matrix?> {
            val n = brs.rank()

            val s = brs.simpleReflections
            val wb = brs.twoColorDD(true)
            val wv = wb.first
            val bv = wb.second

            val w = wv.fold(Matrix(n), {m, i -> m.mul(s[i])})
            val b = bv.fold(Matrix(n), {m, i -> m.mul(s[i])})

            var ww = Matrix(n)
            var wwo: Matrix
            var len: Int
            var newLen = 0
            var flag = true

            do {
                len = newLen
                wwo = ww
                if (flag) {
                    ww = ww.mul(w)
                    flag = false
                } else {
                    ww = ww.mul(b)
                    flag = true
                }

                newLen = brs.length(ww)
            } while (newLen > len)

            var ww2: Matrix? = null
            for (v in 0..(n-1)) {
                val ww2c = wwo.mul(s[v])
                val len2 = brs.length(ww2c)
                if (len2 == len - 1) {
                    ww2 = ww2c
                    break
                }
            }

            System.out.println("WWO: "+brs.myPositiveRoots.intersect(sim.map { wwo.mul(it) }).size)
            System.out.println("WW2: "+brs.myPositiveRoots.intersect(sim.map { ww2?.mul(it) }).size)
            System.out.println("WWO length: "+brs.length(wwo))
            System.out.println("WW2 length: "+brs.length(ww2!!))


            return Pair(wwo, ww2)
        }
    }

    init {
        for (i in 0 until basis.height())
            if (i != noReflection)
                carrier.add(Vector(basis.myCoo[i]).getReflection())
        val iter = 0
        var flag = false
        var newElements: MutableSet<Matrix>
        do {
            newElements = HashSet()

            var count = carrier.size
            //int oldCount = count;
            val c2 = ArrayList(carrier)
            for (j in c2.indices.reversed()) {
                val e = c2[j]
                for (f in carrier) {
                    val g = e.mul(f)
                    if (!carrier.contains(g) && !newElements.contains(g)) {
                        newElements.add(g)
                        count++
                        if (count >= cutoff) {
                            println("Finished @ $count")
                            flag = true
                            break
                        }
                    }
                }
                if (flag) break
            }

            carrier.addAll(newElements)

            //System.out.println("Iteration #" + (iter++) + "; found total # of Weyl group elements: " + carrier.size());
        } while (newElements.size > 0)
        carrier.addAll(newElements)
    }
}
