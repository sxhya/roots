class Matrix {
    val myCoo: Array<DoubleArray>

    constructor(coo: Array<DoubleArray>) {
        myCoo = Array(coo.size) { i -> DoubleArray(coo[i].size)}
        for (i in myCoo.indices) {
            System.arraycopy(coo[i], 0, myCoo[i], 0, myCoo[i].size)
        }
    }

    constructor(m: Matrix): this(m.myCoo)

    constructor(n: Int) { // Identity matrix constructor
        myCoo = Array(n) { DoubleArray(n) }
        for (i in 0 until n) myCoo[i][i] = 1.0
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

    fun column(i: Int): Vector = Vector(DoubleArray(height()) { j -> myCoo[j][i]})

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

    fun mul(m2: Matrix): Matrix {
        if (m2.height() != width()) throw IllegalArgumentException("Wrong dimension of matrices")
        val result = Array(height()) { DoubleArray(m2.width()) }
        for (k in 0 until m2.width()) {
            val c = mul(m2.column(k)).myCoo
            for (i in 0 until height()) result[i][k] = c[i]
        }

        return Matrix(result)
    }

    fun solve (v: Vector): Vector {
        val n: Int = width()
        if (n != height()) throw IllegalArgumentException()
        val am = Matrix(this)
        val a = am.myCoo
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

    fun transpose(): Matrix {
        val w: Int = width()
        val h: Int = height()
        val coo = Array(w) { DoubleArray(h) }
        for (i in 0 until w) for (j in 0 until h) coo[i][j] = myCoo[j][i]
        return Matrix(coo)
    }

    fun inverse(): Matrix {
        val n = this.height()
        assert(n == this.width())

        val result = Array(n) { DoubleArray(n) }
        for (i in 0 until n) {
            val z = DoubleArray(n)
            for (j in 0 until n) z[j] = 0.0
            z[i] = 1.0
            val coli = solve(Vector(z))
            for (j in 0 until n) result[j][i] = coli.myCoo[j]
        }

        return Matrix(result)
    }

    override fun equals(other: Any?): Boolean {
        if (other is Matrix) {
            val m = other as Matrix?
            if (m!!.myCoo.size != myCoo.size) return false
            for (i in myCoo.indices) {
                if (m.myCoo[i].size != myCoo[i].size) return false
                var norm = 0.0
                for (j in 0 until myCoo[i].size) {
                    val diff = myCoo[i][j] - m.myCoo[i][j]
                    norm += diff * diff
                }
                if (Math.sqrt(norm) >= Vector.EPS) return false
            }
            return true
        }
        return false
    }

    override fun hashCode(): Int {
        var result = 0
        for (aMyCoo in myCoo)
            for (anAMyCoo in aMyCoo) {
                val k = Math.round(anAMyCoo / Vector.EPS).toInt()
                result = result * 31 + k
            }
        return result
    }

    override fun toString(): String {
        val result = StringBuilder()
        for (i in myCoo.indices) {
            val str = StringBuilder()
            for (j in 0 until myCoo[i].size) {
                val s = Vector.fmt(myCoo[i][j])
                str.append(s)
                if (j < myCoo[i].size - 1) str.append(", ")
            }
            result.append(str)
            if (i < myCoo.size - 1) result.append("\n")
        }
        return result.toString()
    }

    fun add(m1: Matrix, coo: Double = 1.0): Matrix {
        assert (m1.width() == width())
        assert (m1.height() == height())
        return Matrix(Array(height()) { i -> DoubleArray(width()) { j -> coo*m1.myCoo[i][j] + myCoo[i][j]} })
    }

    companion object {

        fun add(a: Array<DoubleArray>, b: DoubleArray, s: Int, d: Int, xi: Double) {
            val n = a.size
            for (k in 0 until n) a[d][k] += xi * a[s][k]
            b[d] += xi * b[s]
        }

        fun random(size: Int, scale: Double, round: Boolean): Matrix {
            val result = Array(size) { DoubleArray(size) }
            for (i in 0 until size)
                for (j in 0 until size) {
                    result[i][j] = Math.random() * (scale - 0.5)
                    if (round) result[i][j] = Math.round(result[i][j]).toDouble()
                }
            for (i in 0 until size) for (j in 0 until size) if (result[i][j] == 0.0) result[j][i] = 0.0
            return Matrix(result)
        }
    }

}
