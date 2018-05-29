class WeightDiagram(val brs: BasedRootSystem, i: Int) {
    val myWeights: Set<Vector> = brs.weightOrbit(Vector(brs.myFundamentalWeights.myCoo[i]))
    val myEdges: MutableMap<Pair<Vector, Vector>, Int> = HashMap()
    val myWeightsNumbers: MutableMap<Vector, Int> = HashMap()
    val myWeightsOrdered: MutableList<Vector> = ArrayList()

    init {
        //Somehow we should also process zero weights

        for (source in myWeights)
            for (target in myWeights) {
                val d = Vector.minus(source, target)
                for (j in 0 until brs.rank())
                    if (d == Vector(brs.myBasis.myCoo[j]))
                        myEdges[Pair(source, target)] = j
            }


        var counter = 0;
        val edgesCopy = HashMap(myEdges).map { it.key }.toMutableSet()
        val noIncoming = myWeights.minus(edgesCopy.map { it.second }).toMutableSet()

        while (noIncoming.isNotEmpty()) {
            val n = noIncoming.first()
            noIncoming.remove(n)
            myWeightsNumbers[n] = counter
            counter++
            myWeightsOrdered.add(n)
            val edgesToRemove = HashSet<Pair<Vector, Vector>>()
            for (e in edgesCopy) if (e.first == n) edgesToRemove.add(e)
            edgesCopy.removeAll(edgesToRemove)
            for (m in edgesToRemove.map { it.second })
                if (edgesCopy.none { it.second == m })
                    noIncoming.add(m)
        }

        if (!edgesCopy.isEmpty()) throw IllegalArgumentException()

    }

    enum class ComparisonResult{Incomparable, Less, Greater, Equal}

    fun compare(a: Vector, b: Vector): ComparisonResult {
        if (!myWeights.contains(a) || !myWeights.contains(b)) throw IllegalArgumentException()

        if (a == b) return ComparisonResult.Equal
        val d = Vector.minus(b, a)
        val isLess = brs.simpleRootCoefficients(d).map { it >= 0 }.fold(true, {f1, f2 -> f1 && f2}) // a_i <= b_i
        val isGreater = brs.simpleRootCoefficients(d).map { it <= 0 }.fold(true, {f1, f2 -> f1 && f2}) // b_i <= a_i

        if (isLess) return ComparisonResult.Less // b > a
        if (isGreater) return ComparisonResult.Greater // a > b
        return ComparisonResult.Incomparable
    }

    fun leq(a: Vector, b: Vector): Boolean = when (compare(a, b)) {
        ComparisonResult.Less, ComparisonResult.Equal -> true
        else -> false
    }

    fun sup(a: Vector, b: Vector): Vector {
        if (!myWeights.contains(a) || !myWeights.contains(b)) throw IllegalArgumentException()

        val upperBound = HashSet<Vector>()
        for (v in myWeights) {
            if (leq(a, v) && leq(b, v)) upperBound.add(v)
        }


        if (upperBound.isEmpty())
            throw IllegalStateException()

        var flag: Boolean
        do {
            flag = false
            val badBounds = HashSet<Vector>()
            for (u in upperBound) for (v in upperBound)
                    if (u != v && compare(u, v) == ComparisonResult.Less) badBounds.add(v)
            if (badBounds.isNotEmpty()) flag = true
            upperBound.removeAll(badBounds)

        } while (flag)

        if (upperBound.size != 1)
            throw IllegalStateException() // currently does not work for representations with zero roots

        return upperBound.first()
    }

    fun print() {
        for (w in myWeightsOrdered) {
            val targets = myEdges.filter { it.key.first == w }

            for (e in targets) {
                System.out.print(""+myWeightsNumbers[w]+"->"+myWeightsNumbers[e.key.second]+" {"+e.value+"} ")
            }

            System.out.println()
        }
    }

    fun printSpringy() {
        System.out.print("graph.addNodes(")
        myWeightsOrdered.mapIndexed { index, _ ->
            System.out.print("'w$index'")
            if (index < myWeightsOrdered.size -1) System.out.print(", ")}
        System.out.println(")")

        System.out.print("graph.addEdges(")
        myEdges.toList().mapIndexed{ index, it ->
            System.out.print("['w"+myWeightsNumbers[it.first.first]+"', 'w"+myWeightsNumbers[it.first.second]+"', {label: '"+it.second+"'}]")
            if (index < myEdges.size -1) System.out.print(", ")}
        System.out.println(")")
    }

    fun h(lambda: Vector, mu: Vector, depth: Int = 100): Int {
        if (depth == 0)
            throw IllegalStateException()

        if (lambda == mu) return 0
        if (compare(lambda, mu) != ComparisonResult.Incomparable) {
            return brs.rootHeight(Vector.minus(lambda, mu))
        }
        val nu = sup(lambda, mu)
        val result = h(lambda, nu, depth-1) + h(nu, mu, depth-1)
        return result
    }

    private fun phantomProcedure(weight: Vector, root: Vector, fundamentals: List<Vector>, roots: Set<Vector>): Vector {
        if (!myWeights.contains(weight) || !roots.contains(root))
            throw IllegalArgumentException()

        for (i in 0 until brs.myBasis.height()) {
            if (root == fundamentals[i]) {
                val sum = Vector.add(weight, fundamentals[i])
                return if (myWeights.contains(sum)) sum else weight
            }
        }

        for (i in 0 until brs.myBasis.height()) {
            val diff = Vector.minus(root, fundamentals[i])
            if (roots.contains(diff)) return phantomProcedure(phantomProcedure(weight, fundamentals[i], fundamentals, roots), diff, fundamentals, roots)
        }

        throw IllegalArgumentException()

    }

    fun highestWeight(): Vector {
        val highestWeight = myWeightsOrdered[0]
        for (i in 1 until myWeightsOrdered.size)
            if (compare(highestWeight, myWeightsOrdered[i]) != ComparisonResult.Greater)
                throw IllegalStateException()
        return highestWeight
    }

    fun phantom(weight: Vector, root: Vector): Vector {
        if (brs.myPositiveRoots.contains(root)) return phantomProcedure(weight, root, brs.myFundamentalRoots, brs.myPositiveRoots)
        //if (brs.myNegativeRoots.contains(root)) return phantomProcedure(weight, root, brs.myMinusFundamentalRoots, brs.myNegativeRoots)

        throw IllegalArgumentException()
    }

    fun c(lambda: Vector, alpha: Vector): Int {
        if (!myWeights.contains(lambda) || !brs.myRootSet.contains(alpha) || !myWeights.contains(Vector.add(lambda, alpha)))
            throw IllegalArgumentException()

        if (brs.myNegativeRoots.contains(alpha))
            return c(Vector.add(lambda, alpha), Vector.minus(alpha))

        val j = h(phantom(lambda, alpha), lambda) - 1
        return if (j % 2 == 0) 1 else -1
    }

    fun x(alpha: Vector, xi: Double): Matrix {
        val result = Array(myWeights.size, {i -> DoubleArray(myWeights.size, {j -> if (i==j) 1.0 else 0.0 })})
        for (w in myWeights) {
            val wa = Vector.add(w, alpha)
            val wi = myWeightsNumbers[w]!!
            if (myWeights.contains(wa)) {
                val wai = myWeightsNumbers[wa]!!
                result[wai][wi] = xi*c(w, alpha)
            }
        }
        return Matrix(result)
    }



    }