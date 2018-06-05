import Vector.Companion.EPS
import org.junit.Test

import java.util.Arrays
import java.util.HashMap



/**
 * Created by user on 6/15/17.
 */

class Tests {
    @Test
    fun testSNFComputer() {
        val k = arrayOf(longArrayOf(4, 0, 0, 0, 0),
                        longArrayOf(0, 0, 0, 0, 0),
                        longArrayOf(0, 0, 0, 0, 0),
                        longArrayOf(0, 0, 0, 0, 0),
                        longArrayOf(0, 0, 0, 0, 0))

        val im = IntegerMatrix(k)
        val computer = IntegerMatrix.SNFComputer(im)
        computer.doComputeSNF()
        println(Arrays.toString(computer.s.getDiagonalEntries()))
    }

    @Test
    fun testInverse() {
        val m = Matrix.random(4, 5.0, false)
        val mi = m.inverse()
        println(m)
        println()
        println(mi)
        println()
        println(m.mul(mi))
    }


    private fun lemma56(base: Matrix, highestWeight: Int) {
        val brs = BasedRootSystem(base.myCoo)
        val wd = WeightDiagram(brs, highestWeight)
        System.out.println("Number of vertices: "+wd.myWeights.size)
        for (alpha in brs.myRootSet)
            for (lambda in wd.myWeights) // nonzero && lambda + alpha -- also nonzero
                if (wd.myWeights.contains(Vector.add(alpha, lambda)))
                    assert (wd.c(lambda, alpha) == wd.c(Vector.add(lambda, alpha), Vector.minus(alpha)))
        for (alpha in brs.myFundamentalRoots)
            for (beta in brs.myFundamentalRoots)
                for (lambda in wd.myWeights) //TODO: can be zero...
                    if (wd.myWeights.contains(Vector.add(alpha, lambda)) &&
                            wd.myWeights.contains(Vector.add(beta, lambda)) &&
                            wd.myWeights.contains(Vector.add(Vector.add(alpha, beta), lambda)))
                        assert (wd.c(lambda, alpha) * wd.c(Vector.add(lambda, alpha), beta) ==
                        wd.c(lambda, beta) * wd.c(Vector.add(lambda, beta), alpha))
    }

    fun testN(base: Matrix, highestWeight: Int) {
        val brs = BasedRootSystem(base.myCoo)
        val wd = WeightDiagram(brs, highestWeight)

        val xa1 = HashMap<Vector, Matrix>()
        val xam1 = HashMap<Vector, Matrix>()
        val N = HashMap<Pair<Vector, Vector>, Int>()
        for (v in brs.myRootSet) {
            xa1[v] = wd.x(v, 1.0)
            xam1[v] = wd.x(v, -1.0)
            //if (brs.rootHeight(v) == 1) { println(xa1[v]) }
        }


        for (a in brs.myRootSet)
            for (b in brs.myRootSet)
                if (brs.myRootSet.contains(Vector.add(a, b))) {
                    val m1 = xa1[a]!!.mul(xa1[b]!!).mul(xam1[a]!!.mul(xam1[b]!!))
                    when (m1) {
                        xa1[Vector.add(a, b)] -> N[Pair(a, b)] = 1
                        xam1[Vector.add(a, b)] -> N[Pair(a, b)] = -1
                        else -> {
                            println("Can't compute N for " + Arrays.toString(brs.simpleRootCoefficients(a)) + ", " + Arrays.toString(brs.simpleRootCoefficients(b)))
                            throw IllegalStateException() //TODO: Implement computation of N's for non-simply-laced root systems
                        }
                    }
                }

        for (alpha in brs.myRootSet)
            for (beta in brs.myRootSet)
                if (brs.myRootSet.contains(Vector.add(alpha, beta))) {
                    for (gamma in brs.myRootSet)
                        if (brs.myRootSet.contains(Vector.add(beta, gamma)) &&
                                brs.myRootSet.contains(Vector.add(alpha, Vector.add(beta, gamma))))
                        assert (Math.abs(N[Pair(beta, gamma)]!! * N[Pair(alpha, Vector.add(beta, gamma))]!! -
                                N[Pair(Vector.add(alpha, beta), gamma)]!! * N[Pair(alpha, beta)]!!) < EPS) // testSNFComputer formula (6) page 9 -- only for simply-laced root systems



                val gamma = Vector.minus(Vector.add(alpha, beta))
                // testSNFComputer formula (1) page 8
                assert (N[Pair(alpha, beta)]!! == N[Pair(Vector.minus(beta), Vector.minus(alpha))]!!)
                assert (N[Pair(alpha, beta)]!! == - N[Pair(Vector.minus(alpha), Vector.minus(beta))]!!)
                assert (N[Pair(alpha, beta)]!! == - N[Pair(beta, alpha)]!!)

                // testSNFComputer formula (2) page 8
                assert (Math.abs(N[Pair(alpha, beta)]!!.toDouble() / Vector.len2(gamma) - N[Pair(beta, gamma)]!!.toDouble() / Vector.len2(alpha)) < EPS)
                assert (Math.abs(N[Pair(alpha, beta)]!!.toDouble() / Vector.len2(gamma) - N[Pair(gamma, alpha)]!!.toDouble() / Vector.len2(beta)) < EPS)
            }
    }

    @Test
    fun testLemma56() {
        lemma56(RootSystems.albase(3), 0)
        lemma56(RootSystems.albase(4), 0)
        lemma56(RootSystems.dlbase(4), 0)
        lemma56(RootSystems.e6base, 0)
        lemma56(RootSystems.e7base, 6)
    }

    @Test
    fun testStructureConstants() {
        testN(RootSystems.albase(3), 0)
        testN(RootSystems.albase(4), 0)
        testN(RootSystems.dlbase(4), 0)
        testN(RootSystems.e6base, 0)
        testN(RootSystems.e7base, 6)
    }

}
