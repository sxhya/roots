import java.util.ArrayList
import java.util.Collections.singletonList
import java.util.HashSet

class BasedRootSystem internal constructor(basis: Array<DoubleArray>) {
    val myBasis: Matrix = Matrix(basis)
    val myBasisT: Matrix = myBasis.transpose()
    val myFundamentalRoots: List<Vector> = myBasis.myCoo.map { Vector(it) }
    val myMinusFundamentalRoots = myFundamentalRoots.map { Vector.minus(it) }
    val myRootLengths: Array<Double> = myFundamentalRoots.map { Vector.len(it) }.toSet().toTypedArray()

    val myCorootsBasis: Matrix
    val myFundamentalWeights: Matrix
    val myFundamentalWeightsT: Matrix
    val myRootSet: Set<Vector> = Utils.rootSystemByBasis(myBasis)
    val myPositiveRoots: MutableSet<Vector> = HashSet()
    val myNegativeRoots: MutableSet<Vector> = HashSet()
    val simpleReflections: Array<Matrix> = myBasis.myCoo.map { Vector(it).getReflection() }.toTypedArray()

    init {
        for (s in myRootSet) {
            val coo = myBasisT.solve(s)
            var isPositive = false
            for (i in coo.myCoo.indices) {
                if (coo.myCoo[i] > 0) {
                    isPositive = true
                    break
                }
            }
            if (isPositive) myPositiveRoots.add(s) else myNegativeRoots.add(s)
        }

        val coo = myBasis.myCoo.copyOf()
        for (v in 0 until coo.size) {
            val vec = Vector(coo[v])
            coo[v] = Vector.mul(vec, 2/Vector.len2(vec)).myCoo
        }
        myCorootsBasis = Matrix(coo)
        val n = myCorootsBasis.height()
        val fw = Array(n) { _ -> Array(n) { _ -> 0.0} }
        for (v in 0 until coo.size) {
            val rhs = Array(n) { i -> if (i == v) 1.0 else 0.0}
            fw[v] = myCorootsBasis.solve(Vector(rhs.toDoubleArray())).myCoo.toTypedArray()
        }
        myFundamentalWeights = Matrix(fw.map { it.toDoubleArray() }.toTypedArray())
        myFundamentalWeightsT = myFundamentalWeights.transpose()

        if (myRootLengths.size > 2) throw IllegalStateException()
        myRootLengths.sort()
    }

    fun isSimplyLaced(): Boolean = myRootLengths.size == 1

    fun isLong(v: Vector): Boolean = if (myRootSet.contains(v))
        when {
            Math.abs(Vector.len(v) - myRootLengths[0]) < Vector.EPS -> false
            else -> true
        } else throw IllegalArgumentException()


    fun length(weylElement: Matrix): Int = myPositiveRoots.filter { myNegativeRoots.contains(weylElement.mul(it)) }.count()

    fun rank(): Int = myBasis.height()

    fun printDD() = println(identifyDDBonds(myFundamentalRoots))

    fun twoColorDD(firstWhite : Boolean): Pair<Set<Int>, Set<Int>> {
        val white = HashSet<Int>()
        val black = HashSet<Int>()
        val l = myBasis.height()
        if (firstWhite) white.add(0) else black.add(0)
        var hasUncoloredVertices = true
        val bonds = identifyDDBonds(myFundamentalRoots)
        while (hasUncoloredVertices) {
            hasUncoloredVertices = false
            for (i in 0 until l)
                if (!white.contains(i) && !black.contains(i)) {
                    hasUncoloredVertices = true
                    val adjW = adjacentVertices(i, bonds).filter { white.contains(it) }
                    val adjB = adjacentVertices(i, bonds).filter { black.contains(it) }
                    if (adjW.size + adjB.size > 1) throw IllegalStateException()
                    if (adjW.isNotEmpty()) black.add(i)
                    if (adjB.isNotEmpty()) white.add(i)
                }
        }
        return Pair(white, black)
    }

    fun simpleRootCoefficients(root: Vector): Array<Int> = roundArray(myBasisT.solve(root).myCoo.toTypedArray())

    fun fundamentalWeightCoefficents(weight: Vector): Array<Int> = roundArray(myFundamentalWeightsT.solve(weight).myCoo.toTypedArray())

    fun weightOrbit(weight : Vector): Set<Vector> {
        val result = HashSet<Vector>()
        var len: Int
        val s = simpleReflections
        result.add(weight)
        do {
            len = result.size
            val newElements = HashSet<Vector>()
            for (elem in result) for (i in 0 until myBasis.height()) newElements.add(s[i].mul(elem))
            result.addAll(newElements)
        } while (result.size > len)

        return result
    }

    fun rootString(root: Vector): List<Int> {
        if (!myPositiveRoots.contains(root)) throw IllegalArgumentException()

        for (i in 0 until myBasis.height())
            if (root == myFundamentalRoots[i])
                return singletonList(i)

        for (i in 0 until myBasis.height()) {
            val v = Vector.minus(root, myFundamentalRoots[i])
            if (myRootSet.contains(v))
                return singletonList(i) + rootString(v)
        }

        throw IllegalArgumentException()
    }

    fun rootHeight(root: Vector): Int {
        val coo = simpleRootCoefficients(root)
        val allPos = coo.map { it >= 0 }.fold(true, {a, b -> a && b})
        val allNeg = coo.map { it <= 0 }.fold(true, {a, b -> a && b})
        val sum = coo.fold(0) { a, b -> a + b}
        return when {
            allPos -> sum
            allNeg -> -sum
            else -> throw IllegalArgumentException()
        }
    }

    companion object {
        private fun roundArray(a: Array<Double>): Array<Int> =
                Array(a.size) { i ->
                    val d = a[i]
                    assert (Math.abs(d - Math.round(d)) < Vector.EPS)
                    Math.round(d).toInt()
                }

        private fun sqr(d: Double): Double = d * d

        private fun identifyDDBonds(l: List<Vector>): List<Bond> {
            val result = ArrayList<Bond>()
            for (i in l.indices)
                for (j in i + 1 until l.size) {
                    val v1 = Vector(l[i])
                    val v2 = Vector(l[j])
                    val mult = Math.round(sqr(2 * Vector.angle(v1, v2))).toInt()
                    if (mult != 0) result.add(Bond(i, j, mult, Vector.len(v1), Vector.len(v2)))
                }
            return result
        }

        private fun adjacentVertices(i: Int, b: Collection<Bond>): Set<Int> {
            val result = HashSet<Int>()
            for (bn in b) {
                if (bn.source == i) result.add(bn.target)
                if (bn.target == i) result.add(bn.source)
            }
            return result
        }

        private class Bond(val source: Int, val target: Int, val multiplicity: Int = 1, val lenSource: Double, val lenTarget: Double) {
            override fun toString(): String {
                val c = when (multiplicity) {
                    1 -> '-'
                    2 -> '='
                    3 -> '≡'
                    else -> throw IllegalArgumentException()
                }

                val b: Char = when {
                    lenTarget - lenSource > Vector.EPS -> '<'
                    lenSource - lenTarget > Vector.EPS -> '>'
                    else -> '-'
                }

                return("" + (source + 1) + c + b + (target + 1) + " ")
            }
        }
    }


}
