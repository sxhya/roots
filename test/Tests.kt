import A3Subsystems.Companion.calculateSpan
import Tests.Relator.Companion.determineRootPosition
import Vector.Companion.EPS
import Vector.Companion.add
import Vector.Companion.minus
import org.junit.Test

import java.util.Arrays
import java.util.Collections.singletonList
import java.util.HashMap
import kotlin.math.roundToInt


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
        val brs = BasedRootSystem(RootSystems.e7base.myCoo)
        val i = 6
        println(brs.myRootSet.size)
        val sigmaPlus = HashSet<Vector>()
        val sigmaMinus = HashSet<Vector>()
        val delta = HashSet<Vector>()

        val highestWeight = Vector(brs.myFundamentalWeights.myCoo[6])
        val weightOrbit = brs.weightOrbit(highestWeight).toList()
        println(weightOrbit.size)
        for ((index, weight1) in weightOrbit.withIndex()) {
            //if (index > 1) for (weight2 in weightOrbit.subList(0, index - 1)) {
                var coo0 = 0
                var coo1 = 0
                var coom1 = 0
                var coo2 = 0
                var coom2 = 0
            for (root in brs.myRootSet) {
                val prod = Vector.prod(weight1, root)
                if (Math.abs(prod - 0) < EPS) coo0++
                if (Math.abs(prod - 1) < EPS) coo1++
                if (Math.abs(prod - 2) < EPS) coo2++
                if (Math.abs(prod + 1) < EPS) coom1++
                if (Math.abs(prod + 2) < EPS) coom2++
            }
                println("$coom2 $coom1 $coo0 $coo1 $coo2")
        }

        /* for (root in brs.myRootSet) {
            val v = brs.simpleRootCoefficients(root)
            val product = Vector.prod(highestWeight, root)
            println(product)
            when {
                v[i] == 0 ->
                    delta.add(root)
                v[i] == 1 ->
                    sigmaPlus.add(root)
                v[i] == -1 ->
                    sigmaMinus.add(root)
            }
        }

        println(sigmaPlus.size)
        println(delta.size)

        for (sigma in sigmaPlus) {
            var summable = 0
            for (tau in sigmaMinus) {
                if (brs.myRootSet.contains(Vector.add(sigma, tau))) summable++
            }
            println(summable)
        } */
    }

    @Test
    fun testA4() {
        val brs = BasedRootSystem(RootSystems.albase(4).myCoo)
        val rv = brs.myFundamentalRoots[0]
        val fw = Vector(brs.myFundamentalWeights.myCoo[0])
        println(fw)
    }

    @Test
    fun testHypothesis() {
        val brs = BasedRootSystem(RootSystems.e7base.myCoo)
        val i = 6
        val rv = brs.myFundamentalRoots[i]
        val fw = Vector(brs.myFundamentalWeights.myCoo[i])
        val Uk = LinkedHashSet<Vector>()

        for (r in brs.myRootSet) {
            val sc = Vector.prod(r, fw)
            if (Math.abs(sc - 1) < EPS) Uk.add(r)
        }

        for (r2 in brs.myRootSet) {
            for (r3 in Uk) {
                val span = calculateSpan(rv, r2, r3)
                if (span.size == 12) {
                    var sigma = 0
                    for (v in span) if (Math.abs(Vector.prod(v, fw) - 1) < EPS) sigma++
                    print(" $sigma")
                    /* if (sigma == 4) {
                        var canFindAlternativeSubsystem = false
                        for (r4 in brs.myRootSet) {
                            val newSpan = calculateSpan(r2, r3, r4)
                            if (newSpan.size == 12) {
                                sigma = 0
                                for (v in newSpan) if (Math.abs(Vector.prod(v, fw) - 1) < EPS) sigma++
                                if (sigma == 3) {
                                    canFindAlternativeSubsystem = true
                                    break
                                }
                            }
                        }
                        println("Span: ${Vector.prod(r2, r3)} canFindAlternativeSubsystem: $canFindAlternativeSubsystem")
                    } */
                }
            }
        }
    }

    @Test
    fun weightLatticeNF() { /* can determine if weights of a representation span the whole weight lattice */
        val brs = BasedRootSystem(RootSystems.e6base.myCoo)
        val weightDiagram = WeightDiagram(brs, 1)
        val list = ArrayList<LongArray>()
        for (w in weightDiagram.myWeightsOrdered)
            list.add(brs.fundamentalWeightCoefficents(w).map { it.toLong() }.toLongArray())
        val integerMatrix = IntegerMatrix(list.toTypedArray())
        val smithForm = IntegerMatrix.SNFComputer(integerMatrix)
        smithForm.doComputeSNF()
        System.out.println(smithForm.s.toString())
    }

    class Relator(val brs: BasedRootSystem, val a: Vector, val b: Vector, val epsA: Int, val epsB: Int) {
        override fun equals(other: Any?): Boolean {
            if (other is Relator) {
                if (other.a == a && other.b == b && other.epsA == epsA && other.epsB == epsB) return true
                if (other.a == b && other.b == a && other.epsA == epsB && other.epsB == epsA) return true
            }
            return false
        }

        override fun hashCode(): Int {
            val aCode = a.hashCode() * 31 + epsA.hashCode()
            val bCode = b.hashCode() * 31 + epsB.hashCode()
            return aCode xor bCode
        }

        override fun toString(): String {
            val result = Arrays.toString(brs.simpleRootCoefficients(a)) + "=" + this.epsA + " / " +
                    Arrays.toString(brs.simpleRootCoefficients(b)) + "=" + this.epsB
            val position = when (determineRootPosition(brs, a, b)) {
                0 -> "2a + b is a root"
                1 -> "a + 2b is a root"
                2 -> "a + b is a root"
                else -> "a + b is not a root"
            }
            return "$result $position"
        }

        companion object {
            fun determineRootPosition(brs: BasedRootSystem, alpha: Vector, beta: Vector): Int {
                val sum = Vector.add(alpha, beta)
                return if (brs.myRootSet.contains(sum)) when {
                    brs.myRootSet.contains(Vector.add(sum, alpha)) -> 0
                    brs.myRootSet.contains(Vector.add(sum, beta)) -> 1
                    else -> 2
                } else -1
            }
        }

    }

    @Test
    fun commutingWeights() {
        val brs = BasedRootSystem(RootSystems.e6base.myCoo)
        val relators = HashSet<Relator>()
        for (alpha in brs.myRootSet)
            for (beta in brs.myRootSet)
                if (alpha != Vector.minus(beta) && alpha != beta) {
                    for (epsilonAlpha in 0..1)
                        for (epsilonBeta in 0..1) if (epsilonAlpha > 0 || epsilonBeta > 0) {
                            //determine root string
                            val position = determineRootPosition(brs, alpha, beta)
                            val epsSum = epsilonAlpha + epsilonBeta
                            val epsilonCondition = when (position) {
                                0 -> epsilonAlpha + epsSum <= 1
                                1 -> epsilonBeta + epsSum <= 1
                                2 -> epsSum <= 1
                                else -> true

                            }
                            if (epsilonCondition) relators.add(Relator(brs, alpha, beta, epsilonAlpha, epsilonBeta))
                        }
                }

        //for (relator in relators) System.out.println(relator.toString())
        val i = 5
        val fw = brs.myCoweights.myCoo[i]
        val refl = brs.simpleReflections[i]
        val set = HashSet<Vector>()
        set.add(Vector(fw))
        set.add(Vector(refl.mul(Vector(fw))))

        /* val weylGroup = WeylGroup(brs.myBasis, 10000, -1)
        val sset = HashSet<Set<Vector>>()
        for (wg in weylGroup.carrier)
            sset.add(set.map { wg.mul(it) }.toHashSet()) */

        val sset = brs.weightMultiOrbit(set)

        System.out.println("There is ${sset.size} pairs of weights")

        for (relator in relators) {
            System.out.print("Testing...")
            var found = false
            for (commutingPair in sset) { // pair of coweights
                val cPair = commutingPair.toTypedArray()
                assert(cPair.size == 2)
                val cP0 = cPair[0]
                val cP1 = cPair[1]
                val a = Vector.prod(cP0, relator.a) + Vector.prod(cP1, relator.a) + relator.epsA
                val b = Vector.prod(cP0, relator.b) + Vector.prod(cP1, relator.b) + relator.epsB
                val aI = a.roundToInt()
                val bI = b.roundToInt()
                val lossOfPrecision = (Math.abs(aI - a) >= EPS) || (Math.abs(bI - b) >= EPS)
                if (lossOfPrecision)
                    System.err.println("Loss of precision...")
                if (aI <= 0 && bI <= 0 && !lossOfPrecision) {
                    System.out.println("ok; ${relator} ${relator.epsA} -> $aI; ${relator.epsB} -> $bI;")
                    found = true
                    break
                }
            }
            if (!found) {
                System.err.println("failed")
                System.err.println("a:${Arrays.toString(brs.simpleRootCoefficients(relator.a))}, epsA:${relator.epsA} b:${Arrays.toString(brs.simpleRootCoefficients(relator.b))} epsB:${relator.epsB}; relation type: $relator")
                break
            }
        }
    }

    @Test
    fun drawHyperCube() {
        val i = 3
        val brs = BasedRootSystem(RootSystems.clbase(4).myCoo)
        val weightDiagram = WeightDiagram(brs, i)
        val m = HashMap<Vector, Set<Vector>>()
        for (w in weightDiagram.myWeights) {
            val nElems = HashSet<Vector>()
            for (r in brs.myRootSet) if (Vector.prod(r, w) >= 0) nElems.add(r)
            m[w] = nElems
        }

        val fw = brs.myFundamentalWeights.myCoo[i]
        val refl = brs.simpleReflections[i]
        val set = HashSet<Vector>()
        set.add(Vector(fw))
        set.add(Vector(refl.mul(Vector(fw))))

        val weylGroup = WeylGroup(brs.myBasis, 4000, -1)
        val sset = HashSet<Set<Vector>>()
        for (wg in weylGroup.carrier)
            sset.add(set.map { wg.mul(it) }.toHashSet())

        for (orb in sset) {
            for (w in orb) System.out.print("$w ")
            System.out.println()
        }
        System.out.println(sset.size)

        System.out.print("graph.addEdges(")
        sset.toList().mapIndexed { index, it ->
            val itarr = it.toTypedArray()
            System.out.print("['w" + weightDiagram.myWeightsNumbers[itarr[0]] + "', 'w" + weightDiagram.myWeightsNumbers[itarr[1]] + "']")
            if (index < sset.size - 1) System.out.print(", ")
        }
        System.out.println(")")

        System.out.println(Arrays.toString(fw))
        System.out.println(Arrays.toString(refl.mul(Vector(fw)).myCoo))
        weightDiagram.printSpringy()


        for (w in weightDiagram.myWeights)
            for (w2 in weightDiagram.myWeights)
                if (brs.myPositiveRoots.minus(m[w]!!.intersect(m[w2]!!)).size == 1)
                    System.out.println(w.toString() + " " + Arrays.toString(brs.fundamentalWeightCoefficents(w)) + "; " +
                            w2 + " " + Arrays.toString(brs.fundamentalWeightCoefficents(w2)))

    }


    private fun lemma56(base: Matrix, highestWeight: Int) {
        val brs = BasedRootSystem(base.myCoo)
        val wd = WeightDiagram(brs, highestWeight)
        System.out.println("Number of vertices: " + wd.myWeights.size)
        for (alpha in brs.myRootSet)
            for (lambda in wd.myWeights) // nonzero && lambda + alpha -- also nonzero
                if (wd.myWeights.contains(Vector.add(alpha, lambda)))
                    assert(wd.c(lambda, alpha) == wd.c(Vector.add(lambda, alpha), Vector.minus(alpha)))
        for (alpha in brs.myFundamentalRoots)
            for (beta in brs.myFundamentalRoots)
                for (lambda in wd.myWeights) //TODO: can be zero...
                    if (wd.myWeights.contains(Vector.add(alpha, lambda)) &&
                            wd.myWeights.contains(Vector.add(beta, lambda)) &&
                            wd.myWeights.contains(Vector.add(Vector.add(alpha, beta), lambda)))
                        assert(wd.c(lambda, alpha) * wd.c(Vector.add(lambda, alpha), beta) ==
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
                            assert(Math.abs(N[Pair(beta, gamma)]!! * N[Pair(alpha, Vector.add(beta, gamma))]!! -
                                    N[Pair(Vector.add(alpha, beta), gamma)]!! * N[Pair(alpha, beta)]!!) < EPS) // testSNFComputer formula (6) page 9 -- only for simply-laced root systems


                    val gamma = Vector.minus(Vector.add(alpha, beta))
                    // testSNFComputer formula (1) page 8
                    assert(N[Pair(alpha, beta)]!! == N[Pair(Vector.minus(beta), Vector.minus(alpha))]!!)
                    assert(N[Pair(alpha, beta)]!! == -N[Pair(Vector.minus(alpha), Vector.minus(beta))]!!)
                    assert(N[Pair(alpha, beta)]!! == -N[Pair(beta, alpha)]!!)

                    // testSNFComputer formula (2) page 8
                    assert(Math.abs(N[Pair(alpha, beta)]!!.toDouble() / Vector.len2(gamma) - N[Pair(beta, gamma)]!!.toDouble() / Vector.len2(alpha)) < EPS)
                    assert(Math.abs(N[Pair(alpha, beta)]!!.toDouble() / Vector.len2(gamma) - N[Pair(gamma, alpha)]!!.toDouble() / Vector.len2(beta)) < EPS)
                }
    }

    fun isClosed(rs: BasedRootSystem, roots: Set<Vector>): Boolean {
        val list = ArrayList<Vector>(roots)
        var result = true
        for (i in 0 until list.size)
            for (j in i + 1 until list.size) {
                val v = Vector.add(list[i], list[j])
                if (rs.myRootSet.contains(v) && !roots.contains(v)) {
                    val cooi = Arrays.toString(rs.simpleRootCoefficients(list[i]))
                    val cooj = Arrays.toString(rs.simpleRootCoefficients(list[j]))
                    val coov = Arrays.toString(rs.simpleRootCoefficients(v))
                    System.out.println("$cooi $cooj $coov")
                    result = false
                }
            }
        return result
    }

    @Test
    fun testDl() {
        val d5 = BasedRootSystem(RootSystems.dlbase(5).myCoo)
        val root1 = d5.myFundamentalRoots[0]
        val dependentRoots = HashSet<Vector>()
        dependentRoots.add(Vector.minus(root1))
        dependentRoots.add(root1)
        for (r in d5.myRootSet) {
            val v = Vector.add(r, root1)
            if (d5.myRootSet.contains(v)) dependentRoots.add(v)
        }

        val otherRoots = d5.myRootSet.minus(dependentRoots).asSequence().toSet()

        System.out.println(isClosed(d5, otherRoots))
    }

    private fun singleDifferenceExperiment(name: String, brs: BasedRootSystem) {
        for (i in 0 until brs.myFundamentalRoots.size) {
            val ui = HashSet<Vector>()
            for (root in brs.myRootSet) if (brs.simpleRootCoefficients(root)[i] > 0) ui.add(root)
            val experimentResult = Utils.subtractClosure(ui, brs).size == brs.myRootSet.size
            System.out.println("RootSystem: $name i=$i, result=$experimentResult")
            System.out.println()
        }
    }

    @Test
    fun differenceExperiment(){
        singleDifferenceExperiment("A2", BasedRootSystem(RootSystems.albase(2).myCoo))
        singleDifferenceExperiment("A3", BasedRootSystem(RootSystems.albase(3).myCoo))
        singleDifferenceExperiment("A4", BasedRootSystem(RootSystems.albase(4).myCoo))
        singleDifferenceExperiment("A5", BasedRootSystem(RootSystems.albase(5).myCoo))
        singleDifferenceExperiment("A6", BasedRootSystem(RootSystems.albase(6).myCoo))
        singleDifferenceExperiment("A7", BasedRootSystem(RootSystems.albase(7).myCoo))
        singleDifferenceExperiment("D4", BasedRootSystem(RootSystems.dlbase(4).myCoo))
        singleDifferenceExperiment("D5", BasedRootSystem(RootSystems.dlbase(5).myCoo))
        singleDifferenceExperiment("B4", BasedRootSystem(RootSystems.blbase(4).myCoo))
        singleDifferenceExperiment("B5", BasedRootSystem(RootSystems.blbase(5).myCoo))
        singleDifferenceExperiment("C4", BasedRootSystem(RootSystems.clbase(4).myCoo))
        singleDifferenceExperiment("C5", BasedRootSystem(RootSystems.clbase(5).myCoo))
        singleDifferenceExperiment("E6", BasedRootSystem(RootSystems.e6base.myCoo))
        singleDifferenceExperiment("E7", BasedRootSystem(RootSystems.e7base.myCoo))
        singleDifferenceExperiment("E8", BasedRootSystem(RootSystems.e8base.myCoo))
        singleDifferenceExperiment("F4", BasedRootSystem(RootSystems.f4base.myCoo))

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

    @Test
    fun testCurtisTits(){
        val rs = BasedRootSystem(RootSystems.f4base.myCoo)
        fun calculateSpan(a: Vector, b : Vector): MutableSet<Vector> {
            val set = HashSet<Vector>()
            set.add(a); set.add(b)
            Utils.reflectClosure(set)
            return set
        }

        val occupiedRoots = HashSet<Vector>()

        for (a in rs.myFundamentalRoots)
            for (b in rs.myFundamentalRoots)
                occupiedRoots.addAll(calculateSpan(a, b))

        val difference = rs.myRootSet.minus(occupiedRoots)
        val sizes = HashSet<Pair<Int, Set<Boolean>>>()

        for (d in difference) for (c in difference) {
            val s = calculateSpan(c, d).size
            val l = HashSet<Boolean>()
            if (s == 6 || s == 4) {
                l.add(rs.isLong(c))
                l.add(rs.isLong(d))
            }
            sizes.add(Pair(s, l))
        }
        for (s in sizes) print("$s, ")

        println()
        println("CT-presentation roots: ${occupiedRoots.size} out out of ${rs.myRootSet.size}")
    }

    fun ht(root: Vector, brs: BasedRootSystem) = brs.simpleRootCoefficients(root).fold(0, {acc, v -> acc + v})

    fun rootString(brs: BasedRootSystem, root: Vector): List<Int> {
        val result = ArrayList<Int>()
        val coo = brs.simpleRootCoefficients(root)
        val ht = coo.fold(0, {acc, v -> acc + v})

        if (ht <= 0) return rootString(brs, Vector.minus(root))
        if (ht == 1) return singletonList(coo.indexOfFirst { it != 0 } + 1)
        for ((index, fr) in brs.myFundamentalRoots.withIndex()) {
            if (coo[index] == 0) continue
            val diff = minus(root, fr)
            if (brs.myRootSet.contains(diff)) {
                result.add(index + 1)
                result.addAll(rootString(brs, diff))
                break
            }
        }
        return result
    }

    fun nasty(brs: BasedRootSystem, wd: WeightDiagram, weight: Vector, root: Vector): Int {
        val rs = rootString(brs, root)
        var w = weight
        var n = 0
        for (i in rs) {
            val ww = Vector.add(w, brs.myFundamentalRoots[i-1])
            if (wd.myWeights.contains(ww)) {
                w = ww
            } else n++
        }
        return n
    }

    fun c(alpha: Vector, lambda: Vector, brs: BasedRootSystem, wd: WeightDiagram): Int {
        val n = if (brs.myPositiveRoots.contains(alpha)) {
            nasty(brs, wd, lambda, alpha)
        } else {
            nasty(brs, wd, add(lambda, alpha), minus(alpha))
        }

        return if (n % 2 == 0) 1 else -1
    }

    @Test
    fun testRootString() {
        val brs = BasedRootSystem(RootSystems.e7base.myCoo)
        val maxHeight = brs.myPositiveRoots.map { ht(it, brs) }.maxOrNull() ?: 0
        val max = brs.myPositiveRoots.find { ht(it, brs) == maxHeight }!!

        /* for (root in brs.myPositiveRoots) {
            val coo = brs.simpleRootCoefficients(root)
            val ht = Math.abs(coo.fold(0, {acc, v -> acc + v}))
            println(Arrays.toString(coo) + " $ht " + rootString(brs, root))
        } */

        val wd = WeightDiagram(brs, 6)
        val fund = Vector(brs.myFundamentalWeights.myCoo[6])

        /*
        //Test Third Look Proposition 1(1)
        for (w in wd.myWeights) if (wd.myWeights.contains(add(w, max)))
            println(c(max, w, brs, wd))

        println("======")
        //Test Third Look Proposition 1(2)
        for (root in brs.myRootSet) {
            var i : Int = 1

            for (w in wd.myWeights) if (wd.myWeights.contains(add(w, root))) {
                i = i * c(root, w, brs, wd)
            }
            val ht = ht(root, brs) % 2
            val c2 = if (ht == 0) -1 else 1
            println("c: $i; ht: $c2 ${i == c2};")
        } */

        val sigma_1_minus = HashSet<Vector>()
        for (root in brs.myRootSet) {
            val coo = brs.simpleRootCoefficients(root)
            if (coo[6] == -1) sigma_1_minus.add(root)
        }
        println("sigma: ${sigma_1_minus.size}")

        for (w in wd.myWeights) {
            if (w == fund) println("FUND") else {
                val root1 = sigma_1_minus.firstOrNull { Vector.len(minus(add(fund, it), w)) < EPS }
                if (root1 != null) {
                    println(" ${c(root1, fund, brs, wd) == c(minus(root1), add(fund, root1), brs, wd)}")
                } else {
                    var exists = false
                    loop1@for (root1 in sigma_1_minus) for (root2 in sigma_1_minus)
                        if (Vector.len(minus(add(add(fund, root1), root2), w)) < EPS ) {
                            val lambda_1 = add(fund, root1)

                            val c1 = c(root1, fund, brs, wd)
                            val c1m = c(minus(root1), lambda_1, brs, wd)

                            val c2 = c(root2, lambda_1, brs, wd)
                            val c2m = c(minus(root2), w, brs, wd)

                            println("${c1 * c2 == c1m * c2m}")
                            exists = true

                        }
                    if (!exists) {
                        loop1@for (root1 in sigma_1_minus) for (root2 in sigma_1_minus) for (root3 in sigma_1_minus)
                            if (Vector.len(minus(add(add(add(fund, root1), root2), root3), w)) < EPS ) {
                                val lambda_1 = add(fund, root1)
                                val lambda_2 = add(add(fund, root1), root2)

                                if (wd.myWeights.contains(lambda_1) && wd.myWeights.contains(lambda_2)) {
                                    val c1 = c(root1, fund, brs, wd)
                                    val c1m = c(minus(root1), lambda_1, brs, wd)

                                    val c2 = c(root2, lambda_1, brs, wd)
                                    val c2m = c(minus(root2), lambda_2, brs, wd)

                                    val c3 = c(root3, lambda_2, brs, wd)
                                    val c3m = c(minus(root3), w, brs, wd)

                                    println("${c1 * c2 * c3 == c1m * c2m * c3m}")
                                    exists = true
                                } else println("something went wrong")
                            }
                    }
                }
            }
            println("=======")
        }

    }

}
