import java.util.*

class IntegerMatrix(val myCoo: Array<LongArray>) {

    constructor(m: IntegerMatrix): this(Array(m.myCoo.size) { i: Int -> LongArray(m.myCoo[i].size) { j: Int -> m.myCoo[i][j]} })

    constructor(n: Int): this(Array(n) { i: Int -> LongArray(n) { j: Int -> if (i==j) 1 else 0} })

    constructor(n: Int, i: Int, j: Int, xi: Long): this(n) {
        if (i == j || i < 0 || j < 0 || i >= n || j >= n) throw IllegalArgumentException()
        myCoo[i][j] = xi
    }

    fun height(): Int {
        return myCoo.size
    }

    fun width(): Int {
        val w = myCoo[0].size
        for (i in 1 until height())
            if (myCoo[i].size != w)
                throw IllegalStateException("Not uniform width of matrix rows -- not a matrix")
        return w
    }

    fun mul(v: Vector): Vector {
        val w = width()
        if (v.myCoo.size != w) throw IllegalArgumentException("Matrix width should coincide with vector height")
        val result = DoubleArray(height())
        for (i in 0 until height()) {
            var c = 0.0
            for (j in 0 until w) c += myCoo[i][j] * v.myCoo[j]
            result[i] = c
        }
        return Vector(result)
    }

    fun mul(m2: IntegerMatrix): IntegerMatrix {
        val w = width()
        val w2 = m2.width()
        if (m2.height() != w) throw IllegalArgumentException("Wrong dimension of matrices")
        val result = Array(height()) { LongArray(w2) }
        for (k in 0 until w2)
            for (i in 0 until height()) {
                var c: Long = 0
                for (j in 0 until w) c += myCoo[i][j] * m2.myCoo[j][k]
                result[i][k] = c
            }
        return IntegerMatrix(result)
    }

    private fun add(a: Array<DoubleArray>, b: DoubleArray, s: Int, d: Int, xi: Double) {
        val n = a.size
        for (k in 0 until n) a[d][k] += xi * a[s][k]
        b[d] += xi * b[s]
    }

    fun solve(v: Vector): Vector {
        val n = width()
        if (n != height()) throw IllegalArgumentException()
        val am = IntegerMatrix(this)
        val a = Array(n) { DoubleArray(n) }
        for (i in 0 until n)
            for (j in 0 until n)
                a[i][j] = am.myCoo[i][j].toDouble()
        val bv = Vector(v)
        val b = bv.myCoo
        if (n != b.size) throw IllegalArgumentException()
        val result = DoubleArray(n)
        for (i in 0 until n) {
            if (a[i][i] == 0.0) {
                var flag = false
                for (j in i + 1 until n)
                    if (a[j][i] != 0.0) {
                        add(a, b, j, i, 1.0)
                        flag = true
                        break
                    }
                if (!flag) throw IllegalArgumentException()
            }
            for (j in i + 1 until n) add(a, b, i, j, -a[j][i] / a[i][i])
        }
        for (i in n - 1 downTo 0) {
            result[i] = b[i]
            for (j in n - 1 downTo i + 1) result[i] -= result[j] * a[i][j]
            result[i] /= a[i][i]
        }
        return Vector(result)
    }

    fun transpose(): IntegerMatrix {
        val w = width()
        val h = height()
        val coo = Array(w) { LongArray(h) }
        for (i in 0 until w) for (j in 0 until h) coo[i][j] = myCoo[j][i]
        return IntegerMatrix(coo)
    }

    override fun equals(other: Any?): Boolean {
        if (other is IntegerMatrix) {
            val m = other as IntegerMatrix?
            if (m!!.myCoo.size != myCoo.size) return false
            for (i in myCoo.indices) {
                if (m.myCoo[i].size != myCoo[i].size) return false
                for (j in 0 until myCoo[i].size) {
                    if (myCoo[i][j] != m.myCoo[i][j]) return false
                }
            }
            return true
        }
        return false
    }

    override fun hashCode(): Int {
        var result = 0
        for (i in myCoo.indices) {
            for (j in 0 until myCoo[i].size) {
                val k = myCoo[i][j].toInt()
                result = result * 31 + k
            }
        }
        return result
    }

    override fun toString(): String {
        val result = StringBuilder()
        for (i in myCoo.indices) {
            val str = StringBuilder()
            for (j in 0 until myCoo[i].size) {
                val s = Vector.fmt(myCoo[i][j].toDouble())
                str.append(s)
                if (j < myCoo[i].size - 1) str.append(", ")
            }
            result.append(str)
            if (i < myCoo.size - 1) result.append("\n")
        }
        return result.toString()
    }

    fun rowTransform(i: Int, j: Int, coo: Long) {
        if (i < 0 || j < 0 || i >= height() || j >= height())
            throw IllegalArgumentException()
        val w = width()
        for (k in 0 until w) {
            myCoo[i][k] += myCoo[j][k] * coo
        }
    }

    fun swapRows(i: Int, j: Int) {
        if (i < 0 || j < 0 || i >= height() || j >= height())
            throw IllegalArgumentException()
        val w = width()
        for (k in 0 until w) {
            val old = myCoo[j][k]
            myCoo[j][k] = myCoo[i][k]
            myCoo[i][k] = old
        }
    }

    fun colTransform(i: Int, j: Int, coo: Long) {
        if (i < 0 || j < 0 || i >= width() || j >= width())
            throw IllegalArgumentException()
        val h = height()
        for (k in 0 until h) {
            myCoo[k][j] += myCoo[k][i] * coo
        }
    }

    fun swapCols(i: Int, j: Int) {
        if (i < 0 || j < 0 || i >= width() || j >= width())
            throw IllegalArgumentException()
        val h = height()
        for (k in 0 until h) {
            val old = myCoo[k][j]
            myCoo[k][j] = myCoo[k][i]
            myCoo[k][i] = old
        }
    }

    fun getDiagonalEntries(): LongArray {
        val m = Math.min(width(), height())
        val result = LongArray(m)
        for (i in 0 until m) {
            result[i] = myCoo[i][i]
        }
        return result
    }

    class SNFComputer(a: IntegerMatrix) {
        val s = IntegerMatrix(a)
        val aw = a.width()
        val ah = a.height()
        var p = IntegerMatrix(ah)
        var pi = IntegerMatrix(ah)
        var q = IntegerMatrix(aw)
        var qi = IntegerMatrix(aw)
        var isEven: BooleanArray

        private val isInSNF: Boolean
            get() {
                for (k in 0 until ah)
                    for (j in 0 until aw)
                        if (s.myCoo[k][j] != 0L && k != j) return false

                var prev: Long = 1
                for (k in 0 until Math.min(aw, ah)) {
                    val d = s.myCoo[k][k]
                    if (prev == 0L) break
                    if (d % prev != 0L) {
                        if (debugFlag) println("We need to reorder diagonal entries")
                        rowTransform(k - 1, k, 1)
                        return false
                    }
                    prev = s.myCoo[k][k]
                }
                return true
            }

        init {
            isEven = BooleanArray(this.aw)
            for (i in 0 until this.aw) {
                var even = true
                for (j in 0 until this.ah) if (a.myCoo[j][i] % 2 != 0L) even = false
                isEven[i] = even
            }
        }

        fun doComputeSNF() {
            while (!isInSNF) {
                for (i in 0 until Math.min(aw, ah)) {
                    var isInNormalForm: Boolean

                    do {
                        for (k in i + 1 until aw) {
                            var f: Long
                            do f = columnReduce(i, i, k) while (f != 0L)
                        }

                        for (k in i + 1 until ah) {
                            var f: Long
                            do f = rowReduce(i, i, k) while (f != 0L)
                        }

                        isInNormalForm = true

                        if (s.myCoo[i][i] != 0L) {
                            for (k in i + 1 until aw) {
                                val f = s.myCoo[i][k] % s.myCoo[i][i]
                                if (f != 0L) isInNormalForm = false
                            }
                        }

                    } while (!isInNormalForm)

                }
            }

            if (debugFlag) {
                println(s.toString() + '\n')
            }
        }

        private fun colTransform(i: Int, j: Int, coo: Long) { // j <- j + i * coo
            if (debugFlag)
                println("Add columns: " + i + "th *" + coo + " to " + j + "th")

            q.colTransform(i, j, coo)
            s.colTransform(i, j, coo)
            qi.rowTransform(i, j, -coo)

            if (debugFlag) {
                println(s.toString() + '\n')
                println(Arrays.toString(isEven))
            }
        }

        private fun rowTransform(i: Int, j: Int, coo: Long) { // i <- i + j*coo
            if (debugFlag)
                println("Add rows: " + j + "th *" + coo + " to " + i + "th")

            p.rowTransform(i, j, coo)
            s.rowTransform(i, j, coo)
            pi.colTransform(i, j, -coo)
            if (debugFlag) {
                println(s.toString() + '\n')
                println(Arrays.toString(isEven))
            }
        }

        private fun swapCols(i: Int, j: Int) {
            s.swapCols(i, j)
            q.swapCols(i, j)
            qi.swapRows(i, j)
            val f = isEven[i]
            isEven[i] = isEven[j]
            isEven[j] = f

            if (debugFlag) {
                println("Swapped columns #$i and $j")
                println(s.toString() + "\n")
                println(Arrays.toString(isEven))
            }

        }

        private fun swapRows(i: Int, j: Int) {
            s.swapRows(i, j)
            p.swapRows(i, j)
            pi.swapCols(i, j)


            if (debugFlag) {
                println("Swapped rows #$i and $j")
                println(s.toString() + "\n")
            }
        }

        private fun columnReduce(row: Int, col1: Int, col2: Int): Long {
            assert(col1 < col2)

            var x = s.myCoo[row][col1]
            var y = s.myCoo[row][col2]
            val f: Long

            if (Math.abs(x) > Math.abs(y) && y != 0L || x == 0L) {
                swapCols(col1, col2)
            }

            x = s.myCoo[row][col1]
            y = s.myCoo[row][col2]

            if (x == 0L && y == 0L) return 0
            f = y / x

            if (isEven[col2] && !isEven[col1] && f % 2 != 0L) {
                val sign = if (f > 0) 1 else -1
                if (Math.abs(x) == Math.abs(y)) {
                    colTransform(col2, col1, -f)
                } else {
                    colTransform(col1, col2, -sign * (Math.abs(f) + 1))
                }

            } else if (f != 0L) {
                colTransform(col1, col2, -f) // y2 <- y2 - fx1
            }

            return f
        }

        private fun rowReduce(col: Int, row1: Int, row2: Int): Long {
            assert(row1 < row2)

            var x = s.myCoo[row1][col]
            var y = s.myCoo[row2][col]

            if (Math.abs(x) > Math.abs(y) && y != 0L || x == 0L) {
                swapRows(row1, row2)
            }

            x = s.myCoo[row1][col]
            y = s.myCoo[row2][col]

            if (x == 0L && y == 0L) return 0

            val f = y / x

            if (f != 0L) {
                rowTransform(row2, row1, -f) //y2 <- y2 - fx1
            }

            return f
        }

        companion object {
            private val debugFlag = false
        }
    }
}
