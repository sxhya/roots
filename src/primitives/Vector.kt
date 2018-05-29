class Vector constructor (val myCoo: DoubleArray) {
    constructor(len: Int): this(DoubleArray(len))
    constructor(vec: Vector): this(DoubleArray(vec.myCoo.size, {p: Int -> vec.myCoo[p]}))

    fun reflectWRT(y: Vector): Vector = add(y, mul(this, -2 * prod(this, y) / len2(this)))

    override fun equals(o: Any?): Boolean {
        if (o == null) return false
        if (o is Vector) return len(add(this, minus(o))) < EPS
        throw IllegalArgumentException()
    }

    override fun hashCode(): Int = myCoo.fold(0, {result: Int, p: Double -> result * 31 + Math.round(p / EPS).toInt()})

    override fun toString(): String {
        val result = StringBuilder("{")
        for (i in myCoo.indices) {
            result.append(fmt(myCoo[i]))
            if (i < myCoo.size - 1) result.append(", ")
        }
        return result.toString() + "}"
    }

    fun getReflection(): Matrix {
        val n = myCoo.size
        val m = Array(n) { DoubleArray(n) }

        var d = len2(this)
        if (Math.abs(d) < EPS) throw IllegalArgumentException("Can not reflect w.r.t. zero vector.")

        d = -2 / d

        for (i in 0 until n) { //column iterator
            val c = d * myCoo[i]
            for (j in 0 until n) { //row iterator
                m[j][i] = c * myCoo[j]
                if (i == j) m[j][i] += 1.0
            }
        }

        return Matrix(m)
    }


    companion object {
        val EPS = 0.01

        fun fmt(d: Double): String {
            return if (Math.abs(d - d.toLong()) <= 0.00001)
                String.format("%d", d.toLong())
            else
                String.format("%e", d)
        }


        fun add(a: Vector, b: Vector): Vector {
            if (a.myCoo.size != b.myCoo.size) throw IllegalArgumentException()
            val result = Vector(a.myCoo.size)
            for (i in a.myCoo.indices) result.myCoo[i] = a.myCoo[i] + b.myCoo[i]
            return result
        }

        fun mul(a: Vector, c: Double): Vector {
            val result = Vector(a.myCoo.size)
            for (i in a.myCoo.indices) result.myCoo[i] = c * a.myCoo[i]
            return result
        }

        fun minus(a: Vector): Vector = mul(a, -1.0)

        fun minus(a: Vector, b: Vector) = add(a, minus(b))

        fun len2(a: Vector): Double = a.myCoo.fold(0.0, {result: Double, b: Double -> result + b*b})

        fun len(a: Vector): Double = Math.sqrt(len2(a))

        fun prod(a: Vector, b: Vector): Double {
            if (a.myCoo.size != b.myCoo.size) throw IllegalArgumentException()
            return a.myCoo.zip(b.myCoo).fold(0.0, {result: Double, p: Pair<Double, Double> -> result + p.first * p.second }) 
        }

        fun angle(a: Vector, b: Vector): Double = prod(a, b) / (len(a) * len(b))

        fun orthogonal(a: Vector, b: Vector): Boolean = Math.abs(angle(a, b)) < EPS
    }
}